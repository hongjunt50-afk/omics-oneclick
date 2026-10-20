"""一键流水线与单项分析调度：负责落盘、进度回报与取消。"""

from __future__ import annotations

import json
import threading
import time
import traceback
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Callable

import numpy as np
import pandas as pd

from . import plotting
from .data_loader import Dataset
from .deg import differential_expression, pairwise_groups
from .enrichment import deg_genes_for_ora, over_representation
from .genesets import load_library
from .gsea import preranked_gsea, rank_metric
from .pca import run_pca
from .qc import sample_qc, standardize
from .wgcna import run_wgcna

ProgressFn = Callable[[float, str], None]
LogFn = Callable[[str, str], None]


class Cancelled(RuntimeError):
    pass


@dataclass
class AnalysisConfig:
    output_dir: Path
    project_name: str | None = None
    group_a: str | None = None
    group_b: str | None = None
    comparisons: list[tuple[str, str]] = field(default_factory=list)
    padj_cut: float = 0.05
    lfc_cut: float = 1.0
    normalization: str = "quantile"
    ora_libraries: list[str] = field(default_factory=lambda: ["KEGG_2021_Human", "GO_Biological_Process_2023"])
    gsea_library: str = "KEGG_2021_Human"
    gsea_permutations: int = 1000
    run_wgcna: bool = True
    wgcna_max_genes: int = 3000
    wgcna_min_module: int = 30
    wgcna_power: float | None = None
    heatmap_top: int = 60
    seed: int = 20260101


@dataclass
class RunResult:
    project_dir: Path
    figures: dict[str, Path] = field(default_factory=dict)
    tables: dict[str, Path] = field(default_factory=dict)
    log_lines: list[str] = field(default_factory=list)
    summary: dict[str, object] = field(default_factory=dict)
    errors: list[str] = field(default_factory=list)


class Runner:
    """把一次分析的全部产物写进 <output_dir>/OmicsOne_[<项目名>_]<时间戳>/。"""

    def __init__(self, dataset: Dataset, config: AnalysisConfig, *, progress: ProgressFn | None = None, log: LogFn | None = None):
        self.dataset = dataset
        self.config = config
        self._progress = progress or (lambda p, m: None)
        self._log = log or (lambda level, msg: None)
        self._cancel = False
        # 置位表示"允许继续"，清空表示暂停；_check() 在步骤边界等待
        self._gate = threading.Event()
        self._gate.set()

        stamp = datetime.now().strftime("%Y%m%d_%H%M%S")
        slug = _project_slug(config.project_name)
        folder = f"OmicsOne_{slug}_{stamp}" if slug else f"OmicsOne_{stamp}"
        self.project_dir = Path(config.output_dir) / folder
        self.fig_dir = self.project_dir / "figures"
        self.tab_dir = self.project_dir / "tables"
        self.result = RunResult(project_dir=self.project_dir)

    # ---------- 基础 ----------
    def cancel(self) -> None:
        self._cancel = True
        self._gate.set()  # 唤醒可能正在等待的暂停点

    def pause(self) -> None:
        """在当前步骤结束后的下一个检查点停下。"""
        self._gate.clear()

    def resume(self) -> None:
        self._gate.set()

    @property
    def paused(self) -> bool:
        return not self._gate.is_set()

    def _check(self) -> None:
        if self._cancel:
            raise Cancelled("用户已终止分析")
        if not self._gate.is_set():
            self._log("WARNING", "已暂停，等待继续…")
            self._gate.wait()
            if self._cancel:
                raise Cancelled("用户已终止分析")
            self._log("INFO", "已继续")

    def _note(self, level: str, message: str) -> None:
        line = f"[{datetime.now().strftime('%H:%M:%S')}] [{level}] {message}"
        self.result.log_lines.append(line)
        self._log(level, message)

    def _step(self, fraction: float, message: str) -> None:
        self._check()
        self._progress(fraction, message)

    def _write_table(self, name: str, frame: pd.DataFrame) -> Path:
        self.tab_dir.mkdir(parents=True, exist_ok=True)
        path = self.tab_dir / f"{name}.csv"
        frame.to_csv(path, index=True, encoding="utf-8-sig")
        self.result.tables[name] = path
        return path

    def _write_figure(self, name: str, path: Path) -> Path:
        self.result.figures[name] = path
        return path

    # ---------- 主流程 ----------
    def run_all(self, steps: list[str] | None = None) -> RunResult:
        wanted = set(steps or ["qc", "deg", "pca", "enrichment", "gsea", "wgcna"])
        # 富集和 GSEA 都以差异分析结果为输入：单项跑时静默补算 DEG，但不写 DEG 的图表
        need_deg_calc = bool(wanted & {"deg", "enrichment", "gsea"})
        need_deg_output = "deg" in wanted
        self.fig_dir.mkdir(parents=True, exist_ok=True)
        self.tab_dir.mkdir(parents=True, exist_ok=True)

        started = time.time()
        self._note("INFO", f"数据：{self.dataset.summary()}")
        for note in self.dataset.notes:
            self._note("INFO", note)

        expr = self.dataset.expr
        groups = self.dataset.groups

        # 1. 质控 + 标准化
        self._step(0.05, "质控与标准化")
        qc = sample_qc(expr)
        for note in qc.notes:
            self._note("INFO", note)
        self._write_table("01_sample_qc", qc.stats.set_index("sample"))
        self._write_table("02_sample_correlation", qc.correlation)
        self._write_figure("qc_boxplot", plotting.qc_boxplot(expr, self.fig_dir / "01_样本表达分布.png"))
        self._write_figure("qc_heatmap", plotting.qc_correlation_plot(qc.correlation, self.fig_dir / "02_样本相关性.png"))

        expr_norm, norm_note = standardize(expr, self.config.normalization)
        if self.config.normalization != "none":
            self._note("INFO", norm_note)
            self._write_table("03_normalized_expression", expr_norm)
            self._write_figure(
                "qc_boxplot_norm", plotting.qc_boxplot(expr_norm, self.fig_dir / "03_标准化后表达分布.png")
            )
        self.result.summary["qc"] = {
            "samples": expr.shape[1],
            "genes": expr.shape[0],
            "notes": qc.notes,
            "outliers": qc.outliers,
            "normalization": norm_note,
        }

        comparisons = self._resolve_comparisons(groups)
        deg_results: dict[str, pd.DataFrame] = {}

        # 2. 差异分析（富集/GSEA 需要结果表，但没有勾选 DEG 时不落 DEG 图表）
        if need_deg_calc:
            if not need_deg_output:
                self._note("INFO", "富集/GSEA 依赖差异分析结果，已自动补算（不输出 DEG 图表）")
            if groups is None:
                self._note("WARNING", "未识别到分组，跳过差异基因分析")
            else:
                for index, (a, b) in enumerate(comparisons):
                    label = f"{a}_vs_{b}"
                    self._step(0.10 + 0.25 * index / max(len(comparisons), 1), f"差异分析 {label}")
                    try:
                        deg = differential_expression(expr_norm, groups, a, b)
                        deg_results[label] = deg.table
                        sig = deg.significant(padj=self.config.padj_cut, lfc=self.config.lfc_cut)
                        if not need_deg_output:
                            continue
                        self._write_table(f"10_DEG_{label}", deg.table.set_index("gene"))
                        self._write_figure(
                            f"volcano_{label}",
                            plotting.volcano_plot(
                                deg.table,
                                self.fig_dir / f"04_火山图_{label}.png",
                                padj_cut=self.config.padj_cut,
                                lfc_cut=self.config.lfc_cut,
                                title=f"差异基因火山图：{a} vs {b}",
                            ),
                        )
                        self._note(
                            "INFO",
                            f"{label}：显著上调 {int((sig['log2FC'] > 0).sum())} 个，下调 "
                            f"{int((sig['log2FC'] < 0).sum())} 个（padj<{self.config.padj_cut}, |log2FC|≥{self.config.lfc_cut}）",
                        )
                        if not sig.empty:
                            try:
                                top_genes = sig.reindex(sig["pvalue"].sort_values().index)["gene"].tolist()
                                self._write_figure(
                                    f"heatmap_{label}",
                                    plotting.heatmap_plot(
                                        expr_norm,
                                        top_genes,
                                        self.fig_dir / f"05_聚类热图_{label}.png",
                                        groups=groups,
                                        title=f"{label} 差异基因聚类热图（Top {min(len(top_genes), self.config.heatmap_top)}）",
                                        max_genes=self.config.heatmap_top,
                                    ),
                                )
                            except Exception as exc:  # 热图失败不应中断主流程
                                self._note("WARNING", f"热图绘制失败：{exc}")
                    except Exception as exc:
                        self._note("ERROR", f"{label} 差异分析失败：{exc}")
                        self.result.errors.append(f"DEG {label}: {exc}")

            if need_deg_output:
                self.result.summary["deg"] = {
                    label: {
                        "n_up": int(((t["padj"] < self.config.padj_cut) & (t["log2FC"] >= self.config.lfc_cut)).sum()),
                        "n_down": int(((t["padj"] < self.config.padj_cut) & (t["log2FC"] <= -self.config.lfc_cut)).sum()),
                    }
                    for label, t in deg_results.items()
                }

        # 3. PCA
        if "pca" in wanted:
            self._step(0.40, "主成分分析")
            try:
                pca = run_pca(expr_norm, groups=groups)
                for name, frame in pca.to_frames().items():
                    self._write_table(f"20_PCA_{name}", frame)
                self._write_figure(
                    "pca", plotting.pca_plot(pca.scores, pca.explained_ratio, self.fig_dir / "06_主成分分析.png", groups=groups)
                )
                self._note("INFO", f"PCA：PC1 解释 {pca.explained_ratio[0] * 100:.1f}%，PC2 解释 {pca.explained_ratio[1] * 100:.1f}%")
                self.result.summary["pca"] = {
                    "pc1": float(pca.explained_ratio[0]),
                    "pc2": float(pca.explained_ratio[1]),
                }
            except Exception as exc:
                self._note("ERROR", f"PCA 失败：{exc}")
                self.result.errors.append(f"PCA: {exc}")

        # 4. 富集（ORA）
        if "enrichment" in wanted and deg_results:
            for index, (label, table) in enumerate(deg_results.items()):
                self._step(0.50 + 0.10 * index / max(len(deg_results), 1), f"富集分析 {label}")
                for lib_name in self.config.ora_libraries:
                    try:
                        library = load_library(lib_name)
                        genes = deg_genes_for_ora(
                            table, padj=self.config.padj_cut, lfc=self.config.lfc_cut, direction="both"
                        )
                        if len(genes) < 3:
                            self._note("WARNING", f"{label} 显著基因过少（{len(genes)} 个），跳过 {library.label} 富集")
                            continue
                        ora = over_representation(genes, list(expr_norm.index), library)
                        short = _lib_short(lib_name)
                        self._write_table(f"30_ORA_{label}_{short}", ora.table.set_index("term"))
                        self._write_figure(
                            f"ora_bar_{label}_{short}",
                            plotting.enrichment_barplot(
                                ora.table,
                                self.fig_dir / f"07_富集条形图_{label}_{short}.png",
                                title=f"{label} · {library.label} 富集（Top 15）",
                            ),
                        )
                        self._write_figure(
                            f"ora_bubble_{label}_{short}",
                            plotting.enrichment_bubble(
                                ora.table,
                                self.fig_dir / f"08_富集气泡图_{label}_{short}.png",
                                title=f"{label} · {library.label} 富集气泡",
                            ),
                        )
                        top = ora.table.head(3)
                        self._note(
                            "INFO",
                            f"{label} / {library.label}：" + "；".join(f"{r.term}(padj={r.padj:.2g})" for r in top.itertuples()),
                        )
                    except Exception as exc:
                        self._note("WARNING", f"{label} 富集分析（{lib_name}）失败：{exc}")

        # 5. GSEA
        if "gsea" in wanted and deg_results:
            primary_label, primary_table = next(iter(deg_results.items()))
            self._step(0.62, f"GSEA 分析 {primary_label}")
            try:
                library = load_library(self.config.gsea_library)
                metric = rank_metric(primary_table, how="log2FC")
                gsea = preranked_gsea(
                    metric,
                    library,
                    n_permutations=self.config.gsea_permutations,
                    seed=self.config.seed,
                )
                table_out = gsea.table.drop(columns=["leading_edge"], errors="ignore")
                self._write_table(f"40_GSEA_{primary_label}", table_out.set_index("term"))
                self._write_table(
                    "41_GSEA_leading_edge",
                    gsea.table.set_index("term")[["NES", "pvalue", "padj", "leading_edge"]],
                )
                self._write_figure(
                    "gsea_bar",
                    plotting.enrichment_barplot(
                        gsea.table,
                        self.fig_dir / f"09_GSEA条形图_{primary_label}.png",
                        title=f"GSEA · {library.label}（{primary_label}）",
                    ),
                )
                for row in gsea.table.reindex(gsea.table["NES"].abs().sort_values(ascending=False).index).head(3).itertuples():
                    hits = [g for g in str(row.leading_edge).split(",") if g]
                    try:
                        self._write_figure(
                            f"gsea_curve_{row.term}",
                            plotting.gsea_plot(
                                metric,
                                hits,
                                self.fig_dir / f"10_GSEA曲线_{_safe(row.term)}.png",
                                term=row.term,
                                nes=row.NES,
                                pvalue=row.pvalue,
                            ),
                        )
                    except Exception as exc:
                        self._note("WARNING", f"GSEA 曲线绘制失败（{row.term}）：{exc}")
                sig = gsea.significant(padj=0.25)
                self._note("INFO", f"GSEA：{len(gsea.table)} 条通路中 FDR<0.25 的有 {len(sig)} 条")
                self.result.summary["gsea"] = {"terms": int(len(gsea.table)), "n_significant": int(len(sig))}
            except Exception as exc:
                self._note("ERROR", f"GSEA 失败：{exc}")
                self.result.errors.append(f"GSEA: {exc}")

        # 6. WGCNA
        if "wgcna" in wanted and self.config.run_wgcna:
            self._step(0.75, "WGCNA 分析")
            try:
                traits = self._build_traits(groups)
                wgcna = run_wgcna(
                    expr_norm,
                    traits=traits,
                    max_genes=self.config.wgcna_max_genes,
                    beta=self.config.wgcna_power,
                    min_module_size=self.config.wgcna_min_module,
                )
                for note in wgcna.notes:
                    self._note("INFO", note)
                self._write_table("50_WGCNA_power", wgcna.power_table.set_index("power"))
                modules_df = wgcna.modules.rename("module").to_frame()
                self._write_table("51_WGCNA_modules", modules_df)
                self._write_table("52_WGCNA_module_eigengenes", wgcna.module_eigengenes)
                self._write_table("53_WGCNA_module_summary", wgcna.module_summary().set_index("module"))
                if not wgcna.module_trait.empty:
                    self._write_table("54_WGCNA_module_trait", wgcna.module_trait.set_index(["module", "trait"]))

                self._write_figure("wgcna_power", plotting.wgcna_power_plot(wgcna.power_table, self.fig_dir / "11_软阈值选择.png"))
                self._write_figure(
                    "wgcna_tree", plotting.wgcna_dendrogram(wgcna, self.fig_dir / "12_WGCNA聚类树.png")
                )
                if not wgcna.module_trait.empty and wgcna.module_trait["trait"].nunique() > 1:
                    self._write_figure(
                        "wgcna_module_trait",
                        plotting.module_trait_plot(wgcna.module_trait, self.fig_dir / "13_模块性状相关性.png"),
                    )

                hub_frames = []
                for module in wgcna.module_summary()["module"]:
                    hub = wgcna.hub_genes(module, expr_norm, top=30)
                    if not hub.empty:
                        hub.insert(0, "module", module)
                        hub_frames.append(hub)
                if hub_frames:
                    self._write_table("55_WGCNA_hub_genes", pd.concat(hub_frames).set_index(["module", "gene"]))

                self.result.summary["wgcna"] = {
                    "beta": float(wgcna.beta),
                    "modules": int(wgcna.module_summary().shape[0]),
                }
            except Exception as exc:
                self._note("ERROR", f"WGCNA 失败：{exc}\n{traceback.format_exc(limit=2)}")
                self.result.errors.append(f"WGCNA: {exc}")

        elapsed = time.time() - started
        self._step(0.92, "生成报告")
        self.result.summary["elapsed_seconds"] = round(elapsed, 1)
        self._note("INFO", f"全部完成，用时 {elapsed:.1f} 秒")

        log_path = self.project_dir / "run_log.txt"
        log_path.write_text("\n".join(self.result.log_lines), encoding="utf-8")
        (self.project_dir / "summary.json").write_text(
            json.dumps(self.result.summary, ensure_ascii=False, indent=2, default=str), encoding="utf-8"
        )
        self.result.tables["run_log"] = log_path
        self._step(1.0, "完成")
        return self.result

    # ---------- 辅助 ----------
    def _resolve_comparisons(self, groups: pd.Series | None) -> list[tuple[str, str]]:
        if self.config.comparisons:
            return list(self.config.comparisons)
        if self.config.group_a and self.config.group_b:
            return [(self.config.group_a, self.config.group_b)]
        if groups is None:
            return []
        return pairwise_groups(groups)

    def _build_traits(self, groups: pd.Series | None) -> pd.DataFrame:
        if groups is None:
            return pd.DataFrame()
        samples = [str(c) for c in self.dataset.expr.columns]
        labels = groups.reindex(samples).astype(str)
        traits = pd.DataFrame(index=samples)
        for group in pd.unique(labels):
            traits[f"group_{group}"] = (labels == group).astype(float).to_numpy()
        return traits


def _safe(text: str) -> str:
    keep = [c for c in str(text) if c.isalnum() or c in "._-"]
    out = "".join(keep)[:50]
    return out or "term"


# GO 三个分支若不细分会撞名（都取首段就成了 GO），这里给它们各自的后缀
_GO_BRANCH = {"Biological_Process": "BP", "Molecular_Function": "MF", "Cellular_Component": "CC"}


def _lib_short(lib_name: str) -> str:
    """库名转文件名片段：GO 分 BP/MF/CC，其余取首段。"""
    for branch, suffix in _GO_BRANCH.items():
        if branch in lib_name:
            return f"GO_{suffix}"
    return lib_name.split("_")[0]


def _project_slug(name: str | None) -> str:
    """项目名转目录片段：保留中英文与数字，剔除路径分隔符等非法字符。"""
    if not name:
        return ""
    keep = [c for c in str(name) if c.isalnum() or c in "-_"]
    return "".join(keep)[:40]
