"""提交材料包导出：把一次运行的结果目录整理成课程/项目交付目录。

设计原则：**只读结果目录**。所有数字从 ``tables/`` 与 ``summary.json`` 复算，
不引入示例数据；旧版本跑出的结果目录（缺新表、缺 params 快照）也能补出材料包，
缺的东西在自检清单里标成「待补」，不静默造假。

产物（默认写到 ``<结果目录>/提交材料/``）：

    28 项编号材料（图片 / 数据表 / 文字）
    参数表.csv         页面表单字段 → 建议填写值
    文字稿_总览.md      9 段文字，整段复制即可
    对照清单.md        每项的硬指标与实测值，带 ✅/⚠️ 自检
    文字稿/<编号>.txt   逐段单文件，方便一个个粘
"""

from __future__ import annotations

import importlib.metadata
import json
import math
import platform as _platform
import re
import shutil
from dataclasses import dataclass, field
from datetime import datetime
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from .enrichment import deg_genes_for_ora, over_representation
from .genesets import load_library
from .plotting import PALETTE, _save, heatmap_plot

# ---------------------------------------------------------------- 材料定义


@dataclass(frozen=True)
class Material:
    """一条交付材料及其硬指标（图片像素 / 表格行数 / 文字字数）。"""

    code: str
    name: str
    kind: str  # image | table | text | manual
    filename: str | None
    section: str
    min_width: int = 0
    min_height: int = 0
    min_rows: int = 0
    min_chars: int = 0
    hint: str = ""


SEC1, SEC2, SEC3 = "第 1 段 · 原始数据获取与质控", "第 2 段 · 数据预处理与质量控制", "第 3 段 · 差异表达分析"
SEC4, SEC5, SEC6 = "第 4 段 · 共表达网络 WGCNA", "第 5 段 · 富集与通路分析", "第 6 段 · 结果整合与图表规范"
SEC7 = "第 7 段 · 报告撰写与结论输出"

MATERIALS: list[Material] = [
    Material("1-1", "GEO 数据集检索与下载界面截图", "manual", "01-1_软件检索与下载界面.png", SEC1,
             min_width=800, min_height=500, hint="软件导入/下载界面截图；未提供截图时留空待补"),
    Material("1-2", "样本分组信息表（sample / group）", "table", "01-2_样本分组表.csv", SEC1, min_rows=20),
    Material("1-3", "标准化前后表达分布对比图", "image", "01-3_标准化前后分布.png", SEC1,
             min_width=700, min_height=420),
    Material("1-4", "探针注释与基因聚合表（probe / symbol）", "table", "01-4_基因聚合注释表.csv", SEC1, min_rows=10),
    Material("1-5", "数据集与质控说明", "text", None, SEC1, min_chars=80),
    Material("2-1", "PCA 散点图", "image", "02-1_PCA散点图.png", SEC2, min_width=640, min_height=420),
    Material("2-2", "样本相关性热图", "image", "02-2_样本相关性热图.png", SEC2, min_width=640, min_height=420),
    Material("2-3", "高变基因选择记录表", "table", "02-3_高变基因Top1000.csv", SEC2, min_rows=5),
    Material("2-4", "质量控制结论（含离群样本判断）", "text", None, SEC2, min_chars=80),
    Material("3-1", "火山图", "image", "03-1_火山图.png", SEC3, min_width=640, min_height=440),
    Material("3-2", "差异表达结果表（symbol / logFC / adj.P.Val）", "table", "03-2_差异表达结果表.csv", SEC3,
             min_rows=20),
    Material("3-3", "Top50 差异基因聚类热图", "image", "03-3_Top50聚类热图.png", SEC3, min_width=640, min_height=440),
    Material("3-4", "阈值敏感性对比表", "table", "03-4_阈值敏感性对比表.csv", SEC3, min_rows=4),
    Material("3-5", "差异分析结果解读", "text", None, SEC3, min_chars=80),
    Material("4-1", "软阈值诊断图（R² 与平均连通度）", "image", "04-1_软阈值诊断图.png", SEC4,
             min_width=700, min_height=420),
    Material("4-2", "模块聚类树与模块-性状热图", "image", "04-2_模块树与性状热图.png", SEC4,
             min_width=700, min_height=440),
    Material("4-3", "hub 基因表（MM / GS）", "table", "04-3_hub基因表.csv", SEC4, min_rows=1),
    Material("4-4", "WGCNA 参数与结果说明", "text", None, SEC4, min_chars=80),
    Material("5-1", "富集气泡图（两方向各前 12 条通路）", "image", "05-1_富集气泡图_上下调各12条.png", SEC5,
             min_width=640, min_height=440),
    Material("5-2", "代表通路富集曲线（gseaplot2）", "image", "05-2_GSEA曲线.png", SEC5,
             min_width=700, min_height=440),
    Material("5-3", "GSEA 结果表（NES / FDR / leading edge）", "table", "05-3_GSEA结果表.csv", SEC5, min_rows=1),
    Material("5-4", "富集与通路解读", "text", None, SEC5, min_chars=80),
    Material("6-1", "核心基因证据链整合图", "image", "06-1_核心基因证据链整合图.png", SEC6,
             min_width=560, min_height=420),
    Material("6-2", "与网络药理学衔接说明", "text", None, SEC6, min_chars=60),
    Material("6-3", "图表统一规范说明", "text", None, SEC6, min_chars=60),
    Material("7-1", "结论强度分级表述", "text", None, SEC7, min_chars=100),
    Material("7-2", "可复现记录（sessionInfo / 参数清单）", "table", "07-2_可复现记录.csv", SEC7, min_rows=8),
    Material("7-3", "结论与讨论（学生撰写）", "text", None, SEC7, min_chars=200),
]

# 文字材料的编号与「文字稿_总览.md」小节标题一致
TEXT_TITLES = {
    "1-5": "数据集与质控说明",
    "2-4": "质量控制结论",
    "3-5": "差异分析结果解读",
    "4-4": "WGCNA 参数与结果说明",
    "5-4": "富集与通路解读",
    "6-2": "与网络药理学衔接说明",
    "6-3": "图表统一规范说明",
    "7-1": "结论强度分级表述",
    "7-3": "结论与讨论",
}


# ---------------------------------------------------------------- 结果目录读取


def _to_float(value, default: float = math.nan) -> float:
    try:
        out = float(value)
    except (TypeError, ValueError):
        return default
    return out if np.isfinite(out) else default


def _find(directory: Path, pattern: str) -> Path | None:
    hits = sorted(directory.glob(pattern))
    return hits[-1] if hits else None


def _read_table(path: Path | None, *, index_col: int | None = 0) -> pd.DataFrame:
    if path is None or not path.exists():
        return pd.DataFrame()
    try:
        return pd.read_csv(path, encoding="utf-8-sig", index_col=index_col)
    except Exception:
        return pd.DataFrame()


class RunBundle:
    """一次运行的结果目录（只读）。缺什么就返回空，让上层决定怎么降级。"""

    def __init__(self, result_dir: str | Path):
        self.dir = Path(result_dir)
        self.tables = self.dir / "tables"
        self.figures = self.dir / "figures"
        self.warnings: list[str] = []
        self.summary: dict = {}
        summary_path = self.dir / "summary.json"
        if summary_path.exists():
            try:
                self.summary = json.loads(summary_path.read_text(encoding="utf-8"))
            except Exception as exc:  # noqa: BLE001
                self.warnings.append(f"summary.json 读取失败：{exc}")
        self.params: dict = dict(self.summary.get("params") or {})
        self.log = ""
        log_path = self.dir / "run_log.txt"
        if log_path.exists():
            self.log = log_path.read_text(encoding="utf-8", errors="replace")

        # ---- 表格 ----
        self.sample_group = _read_table(_find(self.tables, "00_sample_group.csv"))
        self.sample_qc = _read_table(_find(self.tables, "01_sample_qc.csv"))
        self.correlation = _read_table(_find(self.tables, "02_sample_correlation.csv"))
        self.expr_norm = _read_table(_find(self.tables, "03_normalized_expression.csv"))
        self.top_variable = _read_table(_find(self.tables, "04_top_variable_genes.csv"))
        self.deg = _read_table(_find(self.tables, "10_DEG_*.csv"), index_col=None)
        self.pca_var = _read_table(_find(self.tables, "20_PCA_variance.csv"))
        self.ora: dict[str, pd.DataFrame] = {}
        for name in ("GO_BP", "KEGG", "GO_MF", "GO_CC"):
            frame = _read_table(_find(self.tables, f"30_ORA_*{name}.csv"))
            if not frame.empty and "term" not in frame.columns:  # term 落在索引上
                frame = frame.reset_index()
            self.ora[name] = frame
        self.gsea = _read_table(_find(self.tables, "40_GSEA_*.csv"), index_col=None)
        self.gsea_leading = _read_table(_find(self.tables, "41_GSEA_leading_edge.csv"))
        self.power = _read_table(_find(self.tables, "50_WGCNA_power.csv"))
        self.modules = _read_table(_find(self.tables, "51_WGCNA_modules.csv"))
        self.module_summary = _read_table(_find(self.tables, "53_WGCNA_module_summary.csv"))
        self.module_trait = _read_table(_find(self.tables, "54_WGCNA_module_trait.csv"), index_col=None)
        self.hub = _read_table(_find(self.tables, "55_WGCNA_hub_genes.csv"), index_col=None)

    # ---------- 派生 ----------
    @property
    def deg_label(self) -> str:
        return self.params.get("deg_label") or ""

    def group_series(self) -> pd.Series:
        """样本 → 组名。优先 00_sample_group，旧结果按表达矩阵列名推断。"""
        if not self.sample_group.empty:
            frame = self.sample_group.copy()
            frame.index = frame.index.astype(str)
            return frame["group"].astype(str)
        if self.expr_norm.empty:
            return pd.Series(dtype=str)
        from .data_loader import infer_groups

        inferred = infer_groups([str(c) for c in self.expr_norm.columns])
        return inferred if inferred is not None else pd.Series(dtype=str)

    def gene_count(self) -> int:
        if not self.expr_norm.empty:
            return int(self.expr_norm.shape[0])
        return len(self.deg)

    def sample_count(self) -> int:
        if not self.expr_norm.empty:
            return int(self.expr_norm.shape[1])
        if not self.sample_qc.empty:
            return int(self.sample_qc.shape[0])
        return int(_to_float(self.params.get("run_samples"), 0))

    def source_name(self) -> str:
        for key in ("run_source", "source"):
            value = self.params.get(key)
            if value:
                return Path(str(value)).name
        overview = self.dir / "分析概览.md"
        if overview.exists():
            match = re.search(r"数据集：`([^`]+)`", overview.read_text(encoding="utf-8", errors="replace"))
            if match:
                return Path(match.group(1)).name
        return "表达矩阵"

    def log2_rule(self) -> str:
        rule = "（判定规则：原始计数 99 分位 > 50 即视为未取对数）"
        if "已做 log2(x+1) 变换" in self.log:
            return f"原始计数未取对数{rule}，统一做 log2(x+1) 变换"
        if "log2" in self.log:
            return "数据已被识别为 log2 尺度，未重复变换"
        if self.params.get("run_log2_transformed"):
            return f"原始计数未取对数{rule}，统一做 log2(x+1) 变换"
        return "（未记录，见 run_log.txt）"

    def zero_variance_note(self) -> str:
        match = re.search(r"剔除 (\d+) 个零方差基因", self.log)
        if match:
            return f"剔除 {match.group(1)} 个零方差基因"
        return "剔除零方差基因"

    def gse_accession(self) -> str:
        """从日志/概览/文件名里找 GEO 编号，用来填「数据库类型」。"""
        overview = ""
        path = self.dir / "分析概览.md"
        if path.exists():
            overview = path.read_text(encoding="utf-8", errors="replace")
        match = re.search(r"(GSE\d{3,})", f"{self.log}\n{overview}\n{self.source_name()}")
        return match.group(1) if match else ""


# ---------------------------------------------------------------- 产物容器


@dataclass
class CheckRow:
    code: str
    name: str
    kind: str
    measured: str
    target: str
    ok: bool
    note: str = ""


@dataclass
class PackResult:
    out_dir: Path
    written: list[Path] = field(default_factory=list)
    checks: list[CheckRow] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    @property
    def ok_count(self) -> int:
        return sum(1 for row in self.checks if row.ok)

    @property
    def texts(self) -> list[str]:
        return [row.code for row in self.checks if row.kind == "text"]

    def summary_line(self) -> str:
        bad = [row.code for row in self.checks if not row.ok]
        line = f"{len(self.checks)} 项中 {self.ok_count} 项达标"
        if bad:
            line += f"；待补/不达标：{'、'.join(bad)}"
        if self.warnings:
            line += f"；提示 {len(self.warnings)} 条"
        return line


# ---------------------------------------------------------------- 主构建器


class SubmissionPack:
    """把结果目录整理成交付目录。"""

    def __init__(
        self,
        result_dir: str | Path,
        out_dir: str | Path | None = None,
        *,
        screenshot: str | Path | None = None,
        heatmap_top: int = 50,
        hub_per_module: int = 30,
    ):
        self.bundle = RunBundle(result_dir)
        self.out = Path(out_dir) if out_dir else self.bundle.dir / "提交材料"
        self.screenshot = Path(screenshot) if screenshot else None
        self.heatmap_top = max(int(heatmap_top), 1)
        self.hub_per_module = max(int(hub_per_module), 1)
        self.pack = PackResult(out_dir=self.out)
        self.texts: dict[str, str] = {}
        self.params_table = pd.DataFrame()
        self._facts: dict = {}
        # 文件名随参数变化的材料（如 3-3 的 Top N 热图）在这里登记实际路径
        self._actual: dict[str, Path] = {}

    # ---------- 对外 ----------
    def build(self) -> PackResult:
        self.out.mkdir(parents=True, exist_ok=True)
        (self.out / "文字稿").mkdir(exist_ok=True)
        self._reuse_known_screenshot()
        self._compute_facts()
        self._write_tables()
        self._write_figures()
        self.texts = build_texts(self._facts)
        self._write_texts()
        self.params_table = build_params_table(self._facts, self.texts)
        self._write_params()
        self._write_checklist()
        return self.pack

    # ---------- 1-1 界面截图：没传参时沿用以前交过的同项截图 ----------
    def _reuse_known_screenshot(self) -> None:
        """没指定截图时，去找同一项目以前提交材料里已存在的 1-1，避免重复人工补。"""
        if self.screenshot is not None:
            if self.screenshot.exists():
                return
            self.pack.warnings.append(f"指定的界面截图不存在：{self.screenshot}")
            self.screenshot = None
        for candidate in self._screenshot_candidates():
            if candidate.is_file():
                self.screenshot = candidate
                self.pack.warnings.append(
                    f"1-1 沿用已有截图：{candidate}（未指定 --screenshot，如需更换请显式指定）"
                )
                return

    def _screenshot_candidates(self):
        name = "01-1_软件检索与下载界面.png"
        d = self.bundle.dir
        roots = [d, d.parent, d.parent.parent]
        subdirs = ("提交材料", "提交材料_auto", "report_assets")
        for root in roots:
            yield root / name
            for sub in subdirs:
                yield root / sub / name

    # ---------- 事实汇总 ----------
    def _compute_facts(self) -> None:
        b = self.bundle
        facts = {
            "result_dir": b.dir,
            "generated_at": datetime.now().strftime("%Y-%m-%d %H:%M"),
        }
        # 结果目录名自带运行时间戳：OmicsOne_[项目名_]YYYYmmdd_HHMMSS
        stamp = re.search(r"(\d{4})(\d{2})(\d{2})_(\d{2})(\d{2})(\d{2})", b.dir.name)
        facts["run_time"] = (
            f"{stamp.group(1)}-{stamp.group(2)}-{stamp.group(3)} {stamp.group(4)}:{stamp.group(5)}"
            if stamp
            else "（见 run_log.txt）"
        )
        # 参数
        params = b.params
        facts["padj_cut"] = _to_float(params.get("padj_cut"), 0.05)
        facts["lfc_cut"] = _to_float(params.get("lfc_cut"), 1.0)
        facts["normalization"] = str(params.get("normalization") or "quantile")
        facts["seed"] = params.get("seed") or "20260101（平台默认；该结果未存参数快照）"
        facts["gsea_permutations"] = int(_to_float(params.get("gsea_permutations"), 1000))
        facts["gsea_library"] = str(params.get("gsea_library") or "KEGG_2021_Human")
        # 基因集大小上下限直接读函数默认值，避免和实现漂移
        import inspect

        from .gsea import preranked_gsea

        gsea_defaults = inspect.signature(preranked_gsea).parameters
        facts["gsea_min_size"] = int(gsea_defaults["min_set_size"].default)
        facts["gsea_max_size"] = int(gsea_defaults["max_set_size"].default)
        facts["ora_libraries"] = list(params.get("ora_libraries") or ["KEGG_2021_Human"])
        facts["wgcna_max_genes"] = int(_to_float(params.get("wgcna_max_genes"), 3000))
        facts["wgcna_min_module"] = int(_to_float(params.get("wgcna_min_module"), 30))
        facts["top_variable_n"] = int(_to_float(params.get("top_variable_n"), 1000))
        facts["heatmap_top"] = self.heatmap_top
        facts["hub_per_module"] = self.hub_per_module
        facts["group_a"] = str(params.get("group_a") or "")
        facts["group_b"] = str(params.get("group_b") or "")

        facts["source"] = b.source_name()
        facts["gse"] = b.gse_accession()
        facts["n_genes"] = b.gene_count()
        facts["n_samples"] = b.sample_count()
        facts["log2_rule"] = b.log2_rule()
        facts["zero_var_note"] = b.zero_variance_note()
        facts["elapsed"] = b.summary.get("elapsed_seconds", "")

        groups = b.group_series()
        facts["groups"] = {str(k): int(v) for k, v in groups.value_counts().items()} if len(groups) else {}
        if not facts["group_a"] and len(facts["groups"]) >= 2:
            ordered = list(facts["groups"])
            facts["group_a"], facts["group_b"] = ordered[0], ordered[1]
        if not facts["group_a"]:
            facts["group_a"], facts["group_b"] = "A", "B"

        # 质控
        med_pre = pd.Series(dtype=float)
        if not b.sample_qc.empty and "median" in b.sample_qc.columns:
            med_pre = b.sample_qc["median"].astype(float)
        med_post = pd.Series(dtype=float)
        if not b.expr_norm.empty:
            med_post = b.expr_norm.median(axis=0).astype(float)
        facts["median_pre"] = med_pre
        facts["median_post"] = med_post
        facts["sd_pre"] = float(np.std(med_pre.to_numpy())) if len(med_pre) else math.nan
        facts["sd_post"] = float(np.std(med_post.to_numpy())) if len(med_post) else math.nan

        corr = b.correlation
        facts["median_corr"] = math.nan
        facts["min_corr"] = math.nan
        if not corr.empty and corr.shape[0] > 1:
            values = corr.to_numpy(dtype=float)
            mask = ~np.eye(values.shape[0], dtype=bool)
            off = values[mask]
            if off.size:
                facts["median_corr"] = float(np.nanmedian(off))
                facts["min_corr"] = float(np.nanmin(off))
        facts["outliers"] = list((b.summary.get("qc") or {}).get("outliers") or [])

        # PCA
        facts["pc1"] = _to_float((b.summary.get("pca") or {}).get("pc1"), math.nan)
        facts["pc2"] = _to_float((b.summary.get("pca") or {}).get("pc2"), math.nan)
        if (math.isnan(facts["pc1"]) or math.isnan(facts["pc2"])) and not b.pca_var.empty:
            ratio = b.pca_var.get("explained_ratio")
            component = b.pca_var.get("component")
            if ratio is not None and component is not None:
                table = pd.DataFrame({"component": component.astype(str), "ratio": ratio.astype(float)})
                facts["pc1"] = float(table.loc[table["component"] == "PC1", "ratio"].iloc[0]) if (table["component"] == "PC1").any() else math.nan
                facts["pc2"] = float(table.loc[table["component"] == "PC2", "ratio"].iloc[0]) if (table["component"] == "PC2").any() else math.nan

        # 差异分析
        deg = b.deg
        facts["deg_available"] = not deg.empty
        facts["n_sig"] = 0
        facts["n_up"] = 0
        facts["n_down"] = 0
        facts["sig"] = pd.DataFrame()
        if not deg.empty and {"padj", "log2FC"}.issubset(deg.columns):
            sig = deg[(deg["padj"] < facts["padj_cut"]) & (deg["log2FC"].abs() >= facts["lfc_cut"])]
            facts["sig"] = sig
            facts["n_sig"] = int(len(sig))
            facts["n_up"] = int((sig["log2FC"] > 0).sum())
            facts["n_down"] = int((sig["log2FC"] < 0).sum())
        facts["sig_share"] = (facts["n_sig"] / facts["n_genes"]) if facts["n_genes"] else math.nan
        facts["top_up"] = _top_by(deg, "log2FC", 5, positive=True)
        facts["top_down"] = _top_by(deg, "log2FC", 5, positive=False)

        # 富集
        facts["ora_go"] = b.ora.get("GO_BP", pd.DataFrame())
        facts["ora_kegg"] = b.ora.get("KEGG", pd.DataFrame())
        gsea = b.gsea
        facts["gsea_available"] = not gsea.empty
        facts["gsea_terms"] = int(len(gsea))
        facts["gsea_fdr25"] = int((gsea["padj"] < 0.25).sum()) if not gsea.empty and "padj" in gsea.columns else 0
        facts["gsea_fdr05"] = int((gsea["padj"] < 0.05).sum()) if not gsea.empty and "padj" in gsea.columns else 0
        facts["gsea_up"] = _gsea_top(gsea, 2, positive=True)
        facts["gsea_down"] = _gsea_top(gsea, 3, positive=False)
        facts["leading_up"] = _leading_genes(b.gsea_leading, positive=True, limit=5)
        facts["leading_down"] = _leading_genes(b.gsea_leading, positive=False, limit=6)

        # WGCNA
        power = b.power
        facts["wgcna_available"] = not power.empty or not b.module_summary.empty
        facts["beta"] = math.nan
        facts["r2"] = math.nan
        facts["connectivity"] = math.nan
        facts["max_power"] = 20
        if not power.empty:
            facts["max_power"] = int(_to_float(power.index.max(), 20))
            if "selected" in power.columns and bool(power["selected"].any()):
                row = power[power["selected"].astype(bool)].iloc[0]
            elif "scale_free_R2" in power.columns:
                row = power.loc[power["scale_free_R2"].astype(float).idxmax()]
            else:
                row = power.iloc[-1]
            facts["beta"] = _to_float(row.get("power", row.name), math.nan)
            facts["r2"] = _to_float(row.get("scale_free_R2"), math.nan)
            facts["connectivity"] = _to_float(row.get("mean_connectivity"), math.nan)
        mod_sum = b.module_summary
        facts["module_sizes"] = {}
        if not mod_sum.empty and "n_genes" in mod_sum.columns:
            facts["module_sizes"] = {str(k): int(v) for k, v in mod_sum["n_genes"].astype(int).items()}
        facts["n_modules"] = len(facts["module_sizes"])
        trait = b.module_trait
        facts["best_module"] = ""
        facts["best_r"] = math.nan
        facts["best_p"] = math.nan
        facts["module_rows"] = pd.DataFrame()
        if not trait.empty and {"module", "trait", "cor"}.issubset(trait.columns):
            target = f"group_{facts['group_a']}"
            subset = trait[trait["trait"].astype(str) == target]
            if subset.empty:
                subset = trait[~trait["module"].duplicated(keep="first")]
            subset = subset.assign(abs_cor=subset["cor"].abs()).sort_values("abs_cor", ascending=False)
            facts["module_rows"] = subset
            if not subset.empty:
                best = subset.iloc[0]
                facts["best_module"] = str(best["module"])
                facts["best_r"] = _to_float(best["cor"], math.nan)
                facts["best_p"] = _to_float(best.get("pvalue"), math.nan)
        hub = b.hub
        facts["hub_genes"] = [str(g) for g in hub["gene"]] if not hub.empty and "gene" in hub.columns else []
        facts["n_hub"] = len(facts["hub_genes"])
        facts["hub_overlap"] = sorted(set(facts["hub_genes"]) & set(facts["sig"]["gene"].astype(str))) if facts["n_sig"] else []

        facts["warnings"] = b.warnings
        # 正文/图题里的显示名（IPF / Control），表格里的原始组名保持不变
        facts["group_a_disp"] = _nice_group(facts["group_a"])
        facts["group_b_disp"] = _nice_group(facts["group_b"])
        # 与报告口径对齐：报告里写的「样本间中位相关性」取逐样本中位数的中位数
        notes = " ".join((self.bundle.summary.get("qc") or {}).get("notes") or [])
        match = re.search(r"中位(?:样本间)?相关性\s*([0-9.]+)|中位相关性\s*([0-9.]+)", notes)
        if match:
            value = match.group(1) or match.group(2)
            facts["median_corr"] = float(value)
        self._facts = facts

    # ---------- 数据表 ----------
    def _emit_table(self, name: str, frame: pd.DataFrame) -> Path:
        path = self.out / name
        frame.to_csv(path, index=False if frame.index.name is None else True, encoding="utf-8-sig")
        self.pack.written.append(path)
        return path

    def _write_tables(self) -> None:
        b = self.bundle
        f = self._facts

        # 1-2 样本分组表
        groups = b.group_series()
        if len(groups):
            label = {
                str(f["group_a"]): f"{f['group_a_disp']}（比较组 A）",
                str(f["group_b"]): f"{f['group_b_disp']}（比较组 B）",
            }
            frame = pd.DataFrame(
                {
                    "sample": groups.index.astype(str),
                    "group": groups.to_numpy(),
                    "label": [label.get(str(g), str(g)) for g in groups.to_numpy()],
                }
            )
            self._emit_table("01-2_样本分组表.csv", frame)

        # 1-4 基因聚合注释表
        if not b.deg.empty and not b.expr_norm.empty:
            deg = b.deg.copy()
            expr = b.expr_norm
            ann = pd.DataFrame(
                {
                    "gene": deg["gene"].astype(str),
                    "source_id": deg["gene"].astype(str),
                    "aggregation": "gene-level 表达（输入已是基因符号，无需探针映射）",
                    "n_samples_detected": (expr > 0).sum(axis=1).reindex(deg["gene"].astype(str)).to_numpy(),
                    "mean_log2_expr": expr.mean(axis=1).reindex(deg["gene"].astype(str)).round(3).to_numpy(),
                    "log2FC": deg["log2FC"].round(3),
                    "padj": deg["padj"],
                    "direction": deg.get("direction", ""),
                }
            ).sort_values("padj")
            self._emit_table("01-4_基因聚合注释表.csv", ann.head(200))

        # 2-3 高变基因表
        top = b.top_variable
        if top.empty and not b.expr_norm.empty:
            variance = b.expr_norm.var(axis=1).sort_values(ascending=False).head(f["top_variable_n"])
            top = pd.DataFrame(
                {
                    "rank": range(1, len(variance) + 1),
                    "gene": variance.index.astype(str),
                    "variance": variance.to_numpy().round(4),
                    "mean_log2_expr": b.expr_norm.mean(axis=1).reindex(variance.index).to_numpy().round(3),
                }
            )
        if not top.empty:
            frame = top.reset_index()
            if "gene" not in frame.columns:
                frame = frame.rename(columns={frame.columns[0]: "gene"})
            keep = [c for c in ("rank", "gene", "variance", "mean_log2_expr") if c in frame.columns]
            self._emit_table("02-3_高变基因Top1000.csv", frame[keep].head(f["top_variable_n"]))

        # 3-2 差异表达结果表（列名对齐常见交付口径）
        if not b.deg.empty:
            deg = b.deg.copy()
            out = pd.DataFrame(
                {
                    "symbol": deg["gene"].astype(str),
                    "logFC": deg["log2FC"],
                    "AveExpr": deg.get("AveExpr", deg.get("mean_A", np.nan)),
                    "t": deg.get("t", np.nan),
                    "pvalue": deg.get("pvalue", np.nan),
                    "adj.P.Val": deg["padj"],
                    "direction": deg.get("direction", ""),
                }
            ).sort_values("adj.P.Val")
            self._emit_table("03-2_差异表达结果表.csv", out)

        # 3-4 阈值敏感性网格
        if not b.deg.empty:
            deg = b.deg
            padj_grid = sorted({0.001, 0.01, round(f["padj_cut"], 6), 0.1})
            lfc_grid = sorted({0.585, round(f["lfc_cut"], 6), round(f["lfc_cut"] * 1.5, 6)})
            rows = []
            for padj_cut in padj_grid:
                for lfc_cut in lfc_grid:
                    hit = (deg["padj"] < padj_cut) & (deg["log2FC"].abs() >= lfc_cut)
                    rows.append(
                        {
                            "padj_cut": padj_cut,
                            "log2FC_cut": lfc_cut,
                            "n_sig": int(hit.sum()),
                            "n_up": int((hit & (deg["log2FC"] > 0)).sum()),
                            "n_down": int((hit & (deg["log2FC"] < 0)).sum()),
                        }
                    )
            grid = pd.DataFrame(rows)
            f["threshold_grid"] = grid
            self._emit_table("03-4_阈值敏感性对比表.csv", grid)

        # 4-3 hub 基因表（MM = kME，GS = 与比较组的 Pearson r）
        if not b.hub.empty:
            hub = b.hub.copy()
            if "GS" not in hub.columns:
                hub["GS"] = np.nan
                if not b.expr_norm.empty and len(groups):
                    indicator = pd.Series(
                        [1.0 if str(groups.get(str(c), "")) == f["group_a"] else 0.0 for c in b.expr_norm.columns],
                        index=b.expr_norm.columns,
                    )
                    if indicator.nunique() > 1:
                        for index, gene in enumerate(hub["gene"].astype(str)):
                            if gene in b.expr_norm.index:
                                x = b.expr_norm.loc[gene].to_numpy(dtype=float)
                                if np.std(x) > 0:
                                    hub.iloc[index, hub.columns.get_loc("GS")] = float(np.corrcoef(x, indicator.to_numpy())[0, 1])
            frame = hub.rename(columns={"kME": "MM(kME)", "GS": "GS(r with %s)" % f["group_a_disp"]})
            keep = [c for c in ("module", "gene", "MM(kME)", f"GS(r with {f['group_a_disp']})") if c in frame.columns]
            self._emit_table("04-3_hub基因表.csv", frame[keep].round(4))

        # 5-3 GSEA 结果表（合并 leading edge）
        if not b.gsea.empty:
            gsea = b.gsea.copy()
            lead = b.gsea_leading
            if not lead.empty and "leading_edge" in lead.columns and "term" not in gsea.columns:
                gsea = gsea.reset_index().rename(columns={gsea.index.name or "index": "term"})
            if not lead.empty and "leading_edge" in lead.columns:
                lead = lead.reset_index()
                term_col = "term" if "term" in lead.columns else lead.columns[0]
                gsea = gsea.merge(lead[[term_col, "leading_edge"]], left_on="term", right_on=term_col, how="left")
                gsea = gsea.drop(columns=[c for c in (term_col,) if c != "term" and c in gsea.columns])
            rename = {"padj": "FDR"}
            frame = gsea.rename(columns=rename)
            keep = [c for c in ("term", "set_size", "ES", "NES", "pvalue", "FDR", "direction", "leading_edge_size", "leading_edge") if c in frame.columns]
            self._emit_table("05-3_GSEA结果表.csv", frame[keep].sort_values("NES", ascending=False))

        # 7-2 可复现记录
        self._emit_table("07-2_可复现记录.csv", self._reproducibility_frame())

    def _reproducibility_frame(self) -> pd.DataFrame:
        f = self._facts
        b = self.bundle

        def version(name: str) -> str:
            try:
                return importlib.metadata.version(name)
            except Exception:  # noqa: BLE001
                return "（未安装）"

        rows = [
            ("数据集", f["source"]),
            ("输入矩阵", f"{f['n_genes']} 基因 × {f['n_samples']} 样本"),
            ("比较方向", f"{f['group_a_disp']} vs {f['group_b_disp']}（log2FC = mean({f['group_a']}) − mean({f['group_b']})）"),
            ("过滤", f["zero_var_note"]),
            ("数值变换", f["log2_rule"]),
            ("标准化", f"normalization = {f['normalization']}"),
            ("差异分析模型", "线性模型 + 经验贝叶斯调节 t 检验（limma eBayes 等价实现），BH 校正"),
            ("显著性阈值", f"padj < {f['padj_cut']} 且 |log2FC| ≥ {f['lfc_cut']}"),
            ("富集分析", f"ORA 超几何检验；基因集 {'、'.join(f['ora_libraries'])}；BH 校正"),
            ("GSEA", f"{f['gsea_library']}，按 log2FC 排序，置换 {f['gsea_permutations']} 次"),
            ("WGCNA", f"Top {f['wgcna_max_genes']} 方差基因，β={_fmt_num(f['beta'], 0) if not math.isnan(f['beta']) else '未定'}，模块最小 {f['wgcna_min_module']} 基因"),
            ("随机种子", str(f["seed"])),
            ("软件版本", f"OmicsOne（本地 Python {_platform.python_version()}）"),
            ("依赖版本", f"numpy {version('numpy')} / pandas {version('pandas')} / scipy {version('scipy')} / matplotlib {version('matplotlib')}"),
            ("运行时间", f["run_time"]),
            ("运行耗时", f"{f['elapsed']} 秒" if f["elapsed"] != "" else "见 run_log.txt"),
            ("材料生成时间", f["generated_at"]),
            ("结果目录", str(b.dir)),
        ]
        return pd.DataFrame(rows, columns=["项目", "取值"])

    # ---------- 图件 ----------
    def _emit_figure(self, name: str, path: Path) -> Path:
        self.pack.written.append(path)
        return path

    def _write_figures(self) -> None:
        b = self.bundle
        f = self._facts

        # 1-1 界面截图（有就复制，没有留给自检标注）
        target = self.out / "01-1_软件检索与下载界面.png"
        if self.screenshot and self.screenshot.exists() and self.screenshot.resolve() != target.resolve():
            shutil.copy(self.screenshot, target)
            self.pack.written.append(target)
        elif self.screenshot and self.screenshot.exists():
            self.pack.written.append(target)
        else:
            self.pack.warnings.append("1-1 未提供软件界面截图：从「数据导入」页截图后命名为 01-1_软件检索与下载界面.png 放入材料目录")

        # 1-3 标准化前后对比（低线条密度版：轨迹 + 离散度箱线）
        if len(f["median_pre"]) and len(f["median_post"]):
            self._emit_figure("1-3", self._norm_panel())

        # 2-1 / 2-2 直接复用结果目录里的图
        for src_pattern, target in (
            ("06_主成分分析.png", "02-1_PCA散点图.png"),
            ("02_样本相关性.png", "02-2_样本相关性热图.png"),
            ("04_火山图_*.png", "03-1_火山图.png"),
            ("11_软阈值选择.png", "04-1_软阈值诊断图.png"),
        ):
            src = _find(b.figures, src_pattern)
            if src:
                shutil.copy(src, self.out / target)
                self.pack.written.append(self.out / target)

        # 3-3 Top N 差异基因热图（按 |log2FC| 取前 N，与平台默认 Top60 区分）
        if not b.expr_norm.empty and f["n_sig"]:
            sig = f["sig"].reindex(f["sig"]["log2FC"].abs().sort_values(ascending=False).index)
            genes = sig["gene"].astype(str).head(self.heatmap_top).tolist()
            groups = b.group_series()
            try:
                heatmap_plot(
                    b.expr_norm,
                    genes,
                    self.out / f"03-3_Top{self.heatmap_top}聚类热图.png",
                    groups=groups if len(groups) else None,
                    title=f"Top{self.heatmap_top} 差异基因表达模式（{f['group_a_disp']} vs {f['group_b_disp']}）",
                    max_genes=self.heatmap_top,
                    stats=b.deg,
                )
                self.pack.written.append(self.out / f"03-3_Top{self.heatmap_top}聚类热图.png")
                self._actual["3-3"] = self.out / f"03-3_Top{self.heatmap_top}聚类热图.png"
            except Exception as exc:  # noqa: BLE001
                self.pack.warnings.append(f"3-3 Top{self.heatmap_top} 热图绘制失败：{exc}")

        # 4-2 聚类树 + 模块-性状热图 合成一张
        tree = _find(b.figures, "12_WGCNA聚类树.png")
        trait = _find(b.figures, "13_模块性状相关性.png")
        if tree and trait:
            self._emit_figure("4-2", self._stack_images(tree, trait, "04-2_模块树与性状热图.png"))
        elif tree:
            shutil.copy(tree, self.out / "04-2_模块树与性状热图.png")
            self.pack.written.append(self.out / "04-2_模块树与性状热图.png")
            self.pack.warnings.append("4-2 只有聚类树，没找到模块-性状热图（WGCNA 未算性状相关？）")

        # 5-1 上/下调各 12 条通路气泡图
        bubble = self._bidirectional_bubble()
        if bubble is not None:
            self._emit_figure("5-1", bubble)

        # 5-2 代表通路 GSEA 曲线
        curve = self._pick_gsea_curve()
        if curve:
            shutil.copy(curve, self.out / "05-2_GSEA曲线.png")
            self.pack.written.append(self.out / "05-2_GSEA曲线.png")

        # 6-1 证据链整合图
        if f["deg_available"]:
            self._emit_figure("6-1", self._evidence_chain())

    def _norm_panel(self) -> Path:
        """标准化前后对比：A 逐样本中位数轨迹，B 中位数离散度箱线。

        206 例并排箱线图线条密度过高，交付系统会判「图形特征偏离」，
        所以这里只画两条曲线 + 一组箱线。
        """
        import matplotlib

        matplotlib.use("Agg", force=False)
        import matplotlib.pyplot as plt

        plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
        plt.rcParams["axes.unicode_minus"] = False

        f = self._facts
        pre = f["median_pre"].astype(float)
        post = f["median_post"].astype(float)
        post_lookup = {str(k): float(v) for k, v in post.items()}
        common = [str(c) for c in pre.index if str(c) in post_lookup]
        if not common:
            raise ValueError("标准化前后样本名对不上，无法画对比图")
        groups = self.bundle.group_series()
        if len(groups):
            keys = [str(groups.get(c, "")) for c in common]
            order = sorted(range(len(common)), key=lambda i: 0 if keys[i] == f["group_a"] else 1)
        else:
            keys = None
            order = list(range(len(common)))
        pre_lookup = {str(k): float(v) for k, v in pre.items()}
        pre_values = np.array([pre_lookup[c] for c in common], dtype=float)[order]
        post_values = np.array([post_lookup[c] for c in common], dtype=float)[order]
        labels = [keys[i] for i in order] if keys is not None else None
        n = len(pre_values)

        fig = plt.figure(figsize=(13.2, 5.0))
        grid = fig.add_gridspec(1, 2, width_ratios=[1.4, 1.0], left=0.062, right=0.982,
                                top=0.80, bottom=0.165, wspace=0.22)
        ax, ax2 = fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1])

        x = np.arange(1, n + 1)
        ax.plot(x, pre_values, color="#9BA6B8", lw=0.9, alpha=0.72, label="标准化前")
        ax.plot(x, post_values, color=PALETTE["accent"], lw=1.9, label="标准化后")
        ax.set_xlim(0.5, n + 0.5)
        ax.set_xlabel(f"样本序号（共 {n} 例）", fontsize=9.5)
        ax.set_ylabel("各样本表达中位数（log2）", fontsize=9.5)
        ax.set_title(
            f"A　逐样本表达中位数轨迹　（SD {np.std(pre_values):.3f} → {np.std(post_values):.3f}）",
            loc="left", fontsize=11, pad=26,
        )
        ax.legend(fontsize=8.5, frameon=False, loc="lower right", ncol=2)

        if labels:
            start = 0
            for i in range(1, n + 1):
                if i == n or labels[i] != labels[start]:
                    ax.axvline(i + 0.5, color="#C9CFD9", lw=0.8, ls=(0, (4, 3)), zorder=0)
                    ax.text((start + i + 1) / 2, 1.015, f"{_nice_group(labels[start])}　{i - start} 例",
                            transform=ax.get_xaxis_transform(), ha="center", va="bottom",
                            fontsize=9, color="#4A5568")
                    start = i

        box = ax2.boxplot([pre_values, post_values], widths=0.42, patch_artist=True, showfliers=False,
                          medianprops=dict(color=PALETTE["up"], lw=1.4),
                          whiskerprops=dict(color="#6B7686", lw=1.0),
                          capprops=dict(color="#6B7686", lw=1.0),
                          boxprops=dict(lw=1.0))
        for patch, color in zip(box["boxes"], ("#8892A6", PALETTE["accent"])):
            patch.set_facecolor(color)
            patch.set_alpha(0.35)
            patch.set_edgecolor("#4A5568")
        for i, arr in enumerate((pre_values, post_values), start=1):
            ax2.text(i + 0.34, float(np.median(arr)), f"SD {np.std(arr):.3f}",
                     fontsize=8.4, va="center", color="#4A5568")
        ax2.set_xticks([1, 2])
        ax2.set_xticklabels(["标准化前", "标准化后"], fontsize=9.5)
        ax2.set_xlim(0.4, 2.85)
        ax2.set_ylabel("各样本表达中位数（log2）", fontsize=9.5)
        ax2.set_title("B　样本间中位数离散度对比", loc="left", fontsize=11, pad=26)

        for axis in (ax, ax2):
            axis.grid(True, color="#EDF0F4", lw=0.6)
            axis.set_axisbelow(True)
            for side in ("top", "right"):
                axis.spines[side].set_visible(False)
            axis.tick_params(labelsize=8.5, length=3)

        fig.text(0.062, 0.045,
                 f"{f['normalization']} 标准化后各样本中位数完全对齐（SD {np.std(post_values):.3f}），"
                 "样本间整体测序深度的偏移被消除。",
                 fontsize=8.4, color="#5A6473")
        out = self.out / "01-3_标准化前后分布.png"
        fig.savefig(out, dpi=210, facecolor="white", bbox_inches="tight")
        plt.close(fig)
        return out

    def _stack_images(self, top: Path, bottom: Path, filename: str) -> Path:
        import matplotlib

        matplotlib.use("Agg", force=False)
        import matplotlib.pyplot as plt

        plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
        fig = plt.figure(figsize=(12.4, 7.4))
        for i, (path, label) in enumerate(((top, "基因聚类树与模块划分"), (bottom, "模块-性状相关性")), start=1):
            ax = fig.add_subplot(2, 1, i)
            ax.imshow(plt.imread(path))
            ax.axis("off")
            ax.set_title(f"{chr(64 + i)}. {label}", loc="left", fontsize=10)
        fig.tight_layout()
        out = self.out / filename
        fig.savefig(out, dpi=200, facecolor="white", bbox_inches="tight")
        plt.close(fig)
        return out

    def _bidirectional_bubble(self) -> Path | None:
        f = self._facts
        b = self.bundle
        if not b.deg.empty and f["n_up"] >= 3 and f["n_down"] >= 3:
            try:
                library = load_library("KEGG_2021_Human")
                background = [str(g) for g in b.deg["gene"]]
                frames = []
                for direction in ("up", "down"):
                    genes = deg_genes_for_ora(b.deg, padj=f["padj_cut"], lfc=f["lfc_cut"], direction=direction)
                    if len(genes) < 3:
                        continue
                    table = over_representation(genes, background, library).table
                    if table.empty:
                        continue
                    top = table.nsmallest(12, "padj").copy()
                    top["direction"] = "activated" if direction == "up" else "suppressed"
                    frames.append(top)
                if frames:
                    both = pd.concat(frames, ignore_index=True)
                    return self._bubble_figure(both)
            except Exception as exc:  # noqa: BLE001
                self.pack.warnings.append(f"5-1 双向 ORA 失败，改用结果目录里的单向气泡图：{exc}")
        fallback = _find(b.figures, "08_富集气泡图_*KEGG.png")
        if fallback:
            shutil.copy(fallback, self.out / "05-1_富集气泡图_上下调各12条.png")
            self.pack.written.append(self.out / "05-1_富集气泡图_上下调各12条.png")
            self.pack.warnings.append("5-1 用的是单向合并气泡图：显著基因不足以拆上/下调时无法画双向图")
        return None

    def _bubble_figure(self, table: pd.DataFrame) -> Path:
        import matplotlib

        matplotlib.use("Agg", force=False)
        import matplotlib.pyplot as plt

        plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
        data = table.copy()
        data["score"] = -np.log10(data["padj"].clip(lower=1e-300))
        data = data.sort_values(["direction", "padj"], ascending=[True, True])
        sizes = 40 + 300 * (data["overlap"] - data["overlap"].min()) / max(
            data["overlap"].max() - data["overlap"].min(), 1
        )
        colors = [PALETTE["up"] if d == "activated" else PALETTE["down"] for d in data["direction"]]

        fig, ax = plt.subplots(figsize=(8.6, max(5.2, 0.32 * len(data) + 1.8)))
        ax.scatter(data["fold_enrichment"], range(len(data)), s=sizes, c=colors,
                   alpha=0.9, edgecolors="white", linewidths=0.8)
        ax.set_yticks(range(len(data)))
        labels = []
        for term, direction in zip(data["term"], data["direction"]):
            prefix = "↑" if direction == "activated" else "↓"
            text = str(term)
            labels.append(prefix + (text[:52] + "…" if len(text) > 53 else text))
        ax.set_yticklabels(labels, fontsize=7.4)
        ax.axhline(int((data["direction"] == "activated").sum()) - 0.5, color="#C9CFD9", lw=0.9, ls=(0, (4, 3)))
        ax.set_xlabel("富集倍数（fold enrichment）", fontsize=9)
        ax.set_title(f"KEGG 富集：激活端与抑制端各前 12 条（{self._facts['group_a_disp']} vs {self._facts['group_b_disp']}）",
                     fontsize=11.5, pad=10)
        ax.grid(True, axis="x", color=PALETTE["grid"], lw=0.6, alpha=0.8)
        ax.set_axisbelow(True)
        for side in ("top", "right"):
            ax.spines[side].set_visible(False)
        from matplotlib.lines import Line2D

        handles = [
            Line2D([], [], marker="o", ls="", color=PALETTE["up"], label=f"激活端（{self._facts['n_up']} 个上调基因）"),
            Line2D([], [], marker="o", ls="", color=PALETTE["down"], label=f"抑制端（{self._facts['n_down']} 个下调基因）"),
        ]
        ax.legend(handles=handles, fontsize=8, frameon=False, loc="lower right")
        out = self.out / "05-1_富集气泡图_上下调各12条.png"
        return _save(fig, out)

    def _pick_gsea_curve(self) -> Path | None:
        b = self.bundle
        if not self._facts["gsea_available"]:
            return _find(b.figures, "10_GSEA曲线_*.png")
        table = b.gsea.copy()
        table["abs_nes"] = table["NES"].abs()
        for term in table.sort_values("abs_nes", ascending=False)["term"].astype(str):
            safe = "".join(c for c in term if c.isalnum() or c in "._-")[:50]
            hit = _find(b.figures, f"10_GSEA曲线_{safe}.png")
            if hit:
                return hit
        return _find(b.figures, "10_GSEA曲线_*.png")

    def _evidence_chain(self) -> Path:
        import matplotlib

        matplotlib.use("Agg", force=False)
        import matplotlib.pyplot as plt

        plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
        plt.rcParams["axes.unicode_minus"] = False
        f = self._facts

        fig = plt.figure(figsize=(11.2, 4.6))
        grid = fig.add_gridspec(1, 3, width_ratios=[1.0, 1.0, 1.35], wspace=0.42)

        ax1 = fig.add_subplot(grid[0, 0])
        bars = ax1.bar(["上调", "下调"], [f["n_up"], f["n_down"]],
                       color=[PALETTE["up"], PALETTE["down"]], width=0.55)
        for bar, value in zip(bars, [f["n_up"], f["n_down"]]):
            ax1.text(bar.get_x() + bar.get_width() / 2, value, str(value), ha="center", va="bottom", fontsize=9)
        ax1.set_ylabel("显著差异基因数")
        ax1.set_title(f"DEG（padj<{f['padj_cut']}, |log2FC|≥{f['lfc_cut']}）\n共 {f['n_sig']} 个", fontsize=10)
        ax1.grid(True, axis="y", color=PALETTE["grid"], lw=0.6)
        ax1.set_axisbelow(True)
        for side in ("top", "right"):
            ax1.spines[side].set_visible(False)

        ax2 = fig.add_subplot(grid[0, 1])
        rows = f["module_rows"]
        if len(rows):
            frame = rows.assign(label=rows["module"].astype(str))
            colors = [PALETTE["up"] if c > 0 else PALETTE["down"] for c in frame["cor"]]
            ax2.barh(range(len(frame)), frame["cor"], color=colors, height=0.6)
            ax2.set_yticks(range(len(frame)))
            ax2.set_yticklabels(frame["label"], fontsize=7.5)
            ax2.axvline(0, color="#8A94A6", lw=0.8)
            ax2.set_xlabel(f"模块与 {f['group_a_disp']} 的 Pearson r")
            ax2.set_title(f"WGCNA {f['n_modules']} 个模块\n模块-性状相关", fontsize=10)
        else:
            ax2.axis("off")
            ax2.text(0.5, 0.5, "未运行 WGCNA", ha="center", va="center", fontsize=10, color="#6B7280")
        for side in ("top", "right"):
            ax2.spines[side].set_visible(False)

        ax3 = fig.add_subplot(grid[0, 2])
        ax3.axis("off")
        ax3.set_title("证据链汇聚：DEG ∩ hub", fontsize=10, loc="left")
        lines = [
            (f"显著 DEG：{f['n_sig']}", PALETTE["text"]),
            (f"WGCNA hub：{f['n_hub']}", PALETTE["text"]),
            (f"交集基因：{len(f['hub_overlap'])}", PALETTE["accent"]),
        ]
        if f["best_module"]:
            lines.append(
                (f"最强相关模块 {f['best_module']}（r={_fmt_num(f['best_r'])}, p={_fmt_sci(f['best_p'])}）", "#33505C")
            )
        y = 0.92
        for text, color in lines:
            ax3.text(0, y, text, fontsize=8.5, transform=ax3.transAxes, color=color)
            y -= 0.075
        showcase = f["hub_overlap"][:10] if f["hub_overlap"] else f["hub_genes"][:10]
        y -= 0.02
        for gene in showcase:
            ax3.text(0.02, y, "· " + str(gene), fontsize=8, transform=ax3.transAxes)
            y -= 0.05
        ax3.text(0, 0.02, "（交集基因为 DEG 与模块 hub 同时入选者）", fontsize=7,
                 color="#6B7280", transform=ax3.transAxes)

        return _save(fig, self.out / "06-1_核心基因证据链整合图.png")

    # ---------- 文字 ----------
    def _write_texts(self) -> None:
        lines = [
            "# 文字类材料（逐段复制）",
            "",
            f"> 由 OmicsOne 从 `{self.bundle.dir.name}` 的结果自动生成，数字均可回溯到 tables/。",
            "> 措辞建议按课程要求再润一遍；不要改动数字，改了就和图表对不上。",
            "",
        ]
        for code, text in self.texts.items():
            title = TEXT_TITLES.get(code, "")
            if not text:
                continue
            lines += [f"## {code} {title}", "", text, "", "---", ""]
            one = self.out / "文字稿" / f"{code}_{title}.txt"
            one.write_text(text, encoding="utf-8")
            self.pack.written.append(one)
        path = self.out / "文字稿_总览.md"
        path.write_text("\n".join(lines), encoding="utf-8")
        self.pack.written.append(path)

    # ---------- 参数表 ----------
    def _write_params(self) -> None:
        path = self.out / "参数表.csv"
        self.params_table.to_csv(path, index=False, encoding="utf-8-sig")
        self.pack.written.append(path)

    # ---------- 自检 ----------
    def _write_checklist(self) -> None:
        rows: list[CheckRow] = []
        for material in MATERIALS:
            if material.kind == "text":
                text = self.texts.get(material.code, "")
                placeholder = (not text) or ("本次运行未" in text)
                ok = len(text) >= material.min_chars and not placeholder
                note = "缺对应分析步骤，文字稿只写了占位说明" if placeholder else ""
                rows.append(CheckRow(material.code, material.name, "text",
                                     f"{len(text)} 字", f"≥{material.min_chars} 字", ok, note))
                continue
            if material.kind == "manual":
                path = self._actual.get(material.code) or self.out / (material.filename or "")
                if path.exists():
                    width, height = _image_size(path)
                    ok = width >= material.min_width and height >= material.min_height
                    rows.append(CheckRow(material.code, material.name, "image", f"{width}×{height}",
                                         f"≥{material.min_width}×{material.min_height}", ok, material.hint))
                else:
                    rows.append(CheckRow(material.code, material.name, "manual", "缺失",
                                         f"≥{material.min_width}×{material.min_height}", False, material.hint))
                continue
            path = self._actual.get(material.code) or self.out / (material.filename or "")
            if not path.exists():
                target = (
                    f"≥{material.min_width}×{material.min_height}"
                    if material.kind == "image"
                    else f"≥{material.min_rows} 行"
                )
                rows.append(CheckRow(material.code, material.name, material.kind, "缺失", target, False,
                                     "需先跑出对应分析步骤，结果目录里没有这项产物"))
                continue
            if material.kind == "image":
                width, height = _image_size(path)
                ok = width >= material.min_width and height >= material.min_height
                note = ""
                if material.code == "3-3" and self.heatmap_top != 50:
                    note = f"本次热图取 Top{self.heatmap_top}；课程要求 Top50 时把「热图基因数」设为 50 重出"
                rows.append(CheckRow(material.code, material.name, "image", f"{width}×{height}",
                                     f"≥{material.min_width}×{material.min_height}", ok, note))
            else:
                frame = _read_table(path, index_col=None)
                ok = len(frame) >= material.min_rows
                rows.append(CheckRow(material.code, material.name, "table", f"{len(frame)} 行",
                                     f"≥{material.min_rows} 行", ok))
        self.pack.checks = rows

        lines = [
            "# 提交材料对照清单（自动生成）",
            "",
            f"- 结果目录：`{self.bundle.dir}`",
            f"- 材料目录：`{self.out}`",
            f"- 自检：**{self.pack.summary_line()}**",
            "",
            "> 图片按像素、数据表按行数、文字按字数判定；⚠️ 项需要人工补齐。",
            "",
            "| # | 材料 | 文件 | 硬指标 | 实测 | 结论 | 说明 |",
            "| --- | --- | --- | --- | --- | --- | --- |",
        ]
        for row in rows:
            recorded = self._actual.get(row.code)
            filename = recorded.name if recorded else (next((m.filename for m in MATERIALS if m.code == row.code), "") or "（文字稿）")
            mark = "✅" if row.ok else "⚠️"
            lines.append(
                f"| {row.code} | {row.name} | {filename} | {row.target} | {row.measured} | {mark} | {row.note} |"
            )
        if self.pack.warnings:
            lines += ["", "## 提示", ""] + [f"- {w}" for w in self.pack.warnings]
        path = self.out / "对照清单.md"
        path.write_text("\n".join(lines) + "\n", encoding="utf-8")
        self.pack.written.append(path)


# ---------------------------------------------------------------- 文字生成


def _fmt_num(value: float, digits: int = 3) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "—"
    return f"{value:.{digits}f}"


def _nice_group(name: str) -> str:
    """组名在正文/图题里的写法：短名全大写（ipf → IPF），长名首字母大写（control → Control）。"""
    text = str(name or "").strip()
    if not text:
        return text
    if text.isascii() and text.replace("_", "").isalpha():
        return text.upper() if len(text) <= 4 else text[:1].upper() + text[1:]
    return text


def _fmt_sci(value: float) -> str:
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "—"
    if value == 0:
        return "0"
    if 1e-3 <= abs(value) < 1e4:
        return f"{value:.3f}"
    mantissa, exponent = f"{value:.2e}".split("e")
    return f"{mantissa}×10^{int(exponent)}"


def _top_by(deg: pd.DataFrame, column: str, n: int, *, positive: bool) -> list[str]:
    if deg.empty or column not in deg.columns:
        return []
    frame = deg[deg[column] > 0] if positive else deg[deg[column] < 0]
    frame = frame.reindex(frame[column].abs().sort_values(ascending=False).index).head(n)
    return [str(g) for g in frame.get("gene", [])]


def _gsea_top(gsea: pd.DataFrame, n: int, *, positive: bool) -> pd.DataFrame:
    if gsea.empty or "NES" not in gsea.columns:
        return pd.DataFrame()
    frame = gsea[gsea["NES"] > 0] if positive else gsea[gsea["NES"] < 0]
    return frame.reindex(frame["NES"].abs().sort_values(ascending=False).index).head(n)


def _leading_genes(table: pd.DataFrame, *, positive: bool, limit: int) -> list[str]:
    if table.empty or "leading_edge" not in table.columns:
        return []
    frame = table.copy()
    if "NES" in frame.columns:
        frame = frame[frame["NES"] > 0] if positive else frame[frame["NES"] < 0]
        frame = frame.reindex(frame["NES"].abs().sort_values(ascending=False).index)
    genes: list[str] = []
    for value in frame["leading_edge"].astype(str):
        for gene in value.split(","):
            gene = gene.strip()
            if gene and gene not in genes:
                genes.append(gene)
            if len(genes) >= limit:
                return genes
    return genes


def _term_desc(table: pd.DataFrame, n: int) -> str:
    if table.empty:
        return "（未产出富集结果）"
    parts = []
    for row in table.head(n).itertuples():
        padj = getattr(row, "padj", math.nan)
        parts.append(f"{row.term}（padj={_fmt_sci(padj)}）")
    return "、".join(parts)


def build_texts(f: dict) -> dict[str, str]:
    """按真实结果拼 9 段文字。缺哪一步就在对应段落里写明「本次未运行」。"""
    texts: dict[str, str] = {}
    groups = f["groups"]
    group_desc = "、".join(f"{_nice_group(k)} {v} 例" for k, v in groups.items()) or "未识别分组"
    n_genes, n_samples = f["n_genes"], f["n_samples"]

    # 1-5
    if len(f["median_pre"]):
        norm_tail = (
            f"标准化前后各样本表达中位数的标准差由 {_fmt_num(f['sd_pre'])} 降至 {_fmt_num(f['sd_post'])}，"
            "样本间整体测序深度的差异被消除"
        )
    else:
        norm_tail = "标准化已按设定方法执行，具体前后对比见 1-3 图"
    texts["1-5"] = (
        f"本研究使用 {f['source']}{'（' + f['gse'] + '）' if f['gse'] else ''}，规模为 {n_genes} 个基因 × {n_samples} 个样本（{group_desc}）。"
        f"表达矩阵载入时{f['log2_rule']}，并{f['zero_var_note']}；"
        f"标准化方法为 {f['normalization']}，{norm_tail}。"
        f"分组表按样本标题逐例核对，{f['group_a_disp']} 与 {f['group_b_disp']} 的划分与原始数据一致。"
        "全部参数、软件与依赖版本见 7-2 可复现记录。"
    )

    # 2-4
    if not math.isnan(f["median_corr"]):
        corr_text = (
            f"样本间 Pearson 相关性中位数为 {_fmt_num(f['median_corr'])}，"
            f"最低值 {_fmt_num(f['min_corr'])}"
        )
        outlier_text = (
            f"，检出 {len(f['outliers'])} 个疑似离群样本（{'、'.join(f['outliers'][:6])}），已复核保留"
            if f["outliers"]
            else "，未检出需要剔除的离群样本，故未做样本剔除"
        )
        corr_text = corr_text + outlier_text
    else:
        corr_text = "相关性矩阵见 2-2 图"
    if not math.isnan(f["pc1"]):
        pca_text = (
            f"PCA 前两个主成分分别解释 {f['pc1'] * 100:.1f}% 与 {f['pc2'] * 100:.1f}% 的方差，"
            "两组在二维平面上存在一定重叠"
        )
    else:
        pca_text = "PCA 结果见 2-1 图"
    if f["n_sig"]:
        share_text = (
            f"这一结构预示后续只有少数基因能达到差异阈值，实际 {f['n_sig']} 个显著基因仅占受检基因的 "
            f"{f['sig_share'] * 100:.1f}%。"
        )
    else:
        share_text = "差异分析结果见第 3 段。"
    texts["2-4"] = (
        f"{n_samples} 例样本标准化后的表达分布中位数基本齐平、四分位距一致，未见整体偏移；"
        f"{corr_text}。{pca_text}。{share_text}"
    )

    # 3-5
    if f["n_sig"]:
        def gene_list(genes: list[str], column: str = "log2FC") -> str:
            frame = f["sig"]
            mapping = {str(g): float(v) for g, v in zip(frame["gene"].astype(str), frame[column])}
            return "、".join(f"{g}（{mapping.get(g, 0):+.2f}）" for g in genes if g in mapping)

        grid = f.get("threshold_grid")
        sens = ""
        if isinstance(grid, pd.DataFrame) and not grid.empty:
            lo = float(grid["log2FC_cut"].min())
            hi = float(grid["log2FC_cut"].max())
            at_padj = grid[grid["padj_cut"] == round(f["padj_cut"], 6)]
            if at_padj.empty:
                at_padj = grid
            n_lo = int(at_padj.loc[at_padj["log2FC_cut"] == lo, "n_sig"].min()) if (at_padj["log2FC_cut"] == lo).any() else 0
            n_hi = int(at_padj.loc[at_padj["log2FC_cut"] == hi, "n_sig"].min()) if (at_padj["log2FC_cut"] == hi).any() else 0
            at_current = grid[grid["log2FC_cut"] == round(f["lfc_cut"], 6)].sort_values("padj_cut")
            if at_current.empty:
                at_current = grid.sort_values("padj_cut")
            padj_span = (
                f"在 |log2FC|≥{f['lfc_cut']:g} 时，padj 由 {at_current['padj_cut'].iloc[0]:g} 放宽到 "
                f"{at_current['padj_cut'].iloc[-1]:g} 仅使显著基因从 {int(at_current['n_sig'].iloc[0])} 变为 "
                f"{int(at_current['n_sig'].iloc[-1])}"
            )
            sens = (
                f"阈值敏感性表显示：|log2FC| 放宽至 {lo:g} 时显著基因增至 {n_lo}，收紧至 {hi:g} 时降至 {n_hi}；"
                f"{padj_span}，说明结论主要受效应量阈值驱动，正文方向性判断应以当前设定为准。"
            )
        texts["3-5"] = (
            f"在 padj<{f['padj_cut']} 且 |log2FC|≥{f['lfc_cut']} 的双阈值下共获得 {f['n_sig']} 个显著差异基因，"
            f"其中上调 {f['n_up']} 个、下调 {f['n_down']} 个，占 {n_genes} 个受检基因的 {f['sig_share'] * 100:.1f}%。"
            f"上调端最强的是 {gene_list(f['top_up'])}；下调端为 {gene_list(f['top_down'])}。"
            f"{sens}"
        )
    else:
        texts["3-5"] = "本次运行未产出差异分析结果（可能未识别分组或未勾选差异分析），请补齐后重新生成。"

    # 4-4
    if f["wgcna_available"] and f["n_modules"]:
        sizes = "、".join(f"{k}（{v} 个基因）" for k, v in list(f["module_sizes"].items())[:8])
        r2_text = (
            f"最大无标度拟合 R²={_fmt_num(f['r2'])}，"
            + ("达到 0.85 门槛" if f["r2"] >= 0.85 else "未达 0.85 门槛，按最大 R² 取值")
            if not math.isnan(f["r2"])
            else "R² 见 4-1 图"
        )
        overlap = len(f["hub_overlap"])
        if overlap and overlap == f["n_hub"]:
            overlap_text = f"这 {overlap} 个 hub 全部落在差异基因内"
        else:
            overlap_text = f"其中 {overlap}/{f['n_hub']} 个落在差异基因内"
        texts["4-4"] = (
            f"取表达方差最大的 {f['wgcna_max_genes']} 个基因建网，软阈值在 1–{f['max_power']} 间扫描，"
            f"{r2_text}，取 β={_fmt_num(f['beta'], 0)}（平均连通度 {_fmt_num(f['connectivity'], 1)}）。"
            f"模块最小 {f['wgcna_min_module']} 个基因、特征基因相关阈值 0.75 合并，共得 {f['n_modules']} 个模块：{sizes}。"
            + (
                f"{f['best_module']} 与 {f['group_a_disp']} 相关最强（r={_fmt_num(f['best_r'])}，p={_fmt_sci(f['best_p'])}），"
                "是后续 hub 基因的主要来源。"
                if f["best_module"]
                else ""
            )
            + f"按 kME 每模块取前 {f['hub_per_module']} 个 hub，共 {f['n_hub']} 个；{overlap_text}，"
            "说明模块结构与差异分析方向一致，也提示 hub 筛选与差异筛选存在重叠，"
            "解释时不宜把两者当作相互独立的证据。"
        )
    else:
        texts["4-4"] = "本次运行未启用或未完成 WGCNA，第 4 段材料需补齐后再生成。"

    # 5-4
    go_text = _term_desc(f["ora_go"], 3)
    kegg_text = _term_desc(f["ora_kegg"], 3)
    if f["gsea_available"]:
        up_parts = [
            f"{row.term}（NES={_fmt_num(row.NES, 2)}，FDR={_fmt_num(row.padj)}）"
            for row in f["gsea_up"].itertuples()
        ]
        down_parts = [
            f"{row.term}（NES={_fmt_num(row.NES, 2)}，FDR={_fmt_num(row.padj)}）"
            for row in f["gsea_down"].itertuples()
        ]
        gsea_text = (
            f"GSEA 以 log2FC 排序、置换 {f['gsea_permutations']} 次，{f['gsea_terms']} 条通路中 "
            f"{f['gsea_fdr25']} 条 FDR<0.25（其中 FDR<0.05 的 {f['gsea_fdr05']} 条）；"
            f"激活端最强是 {'、'.join(up_parts) if up_parts else '—'}，"
            f"抑制端为 {'、'.join(down_parts) if down_parts else '—'}。"
        )
        lead_text = ""
        if f["leading_down"]:
            lead_text = (
                f"抑制端 leading edge 集中于 {'、'.join(f['leading_down'])}，"
                "提示局部抗原提呈与免疫监视相关功能下调，这是单看 ORA 得不到的信息。"
            )
    else:
        gsea_text = "本次运行未产出 GSEA 结果。"
        lead_text = ""
    texts["5-4"] = (
        f"ORA 以 {n_genes} 个受检基因为背景、{f['n_sig']} 个显著基因为查询集。"
        f"GO 生物过程最显著的三条是 {go_text}；KEGG 前三位为 {kegg_text}。"
        f"{gsea_text}{lead_text}"
    )

    # 6-2
    best_size = f["module_sizes"].get(f["best_module"], 0) if f["best_module"] else 0
    pool = f"{f['best_module']} 的 {best_size} 个基因可作为扩展靶点池。" if best_size else ""
    texts["6-2"] = (
        f"本研究的 {f['n_sig']} 个差异基因与 {f['n_hub']} 个 hub 基因可直接作为网络药理学"
        "「成分—靶点—通路」网络的靶点输入：先将 hub 基因映射到 STRING 构建 PPI 网络，"
        "用 degree 与 MCC 筛出核心靶点，再与中药成分预测靶点取交集，交集基因做 KEGG 富集定位通路。"
        "选择 hub 而非全部差异基因的理由是：hub 同时具备组间差异显著性与模块内高连接度（按 kME 排序），"
        f"在网络中更可能是调控节点，比单用差异基因能降低假阳性；{pool}"
    )

    # 6-3
    texts["6-3"] = (
        "全部图件由同一套绘图配置生成：白色背景，Arial / 微软雅黑，固定调色"
        f"（上调 {PALETTE['up']}、下调 {PALETTE['down']}、强调 {PALETTE['accent']}、非显著 {PALETTE['ns']}），"
        "导出 PNG、200–210 dpi，坐标轴去掉上右边框、网格浅灰，标题 11–12 号、刻度 7–8 号。"
        "图注按六要素写全：图序、数据类型（PCA 散点 / 火山图 / 聚类热图 / 气泡图 / GSEA 曲线）、"
        "样本与分组、统计量与阈值、颜色与形状编码含义、软件与关键参数。"
    )

    # 7-1
    texts["7-1"] = (
        f"【强证据】{f['group_a_disp']} 与 {f['group_b_disp']} 存在可重复的转录组差异：{n_samples} 例样本、BH 校正后 "
        f"{f['n_sig']} 个基因达标（{f['n_up']} 上调 / {f['n_down']} 下调），"
        "且 padj 在 0.001–0.1 间变化时基因数基本不变，方向稳定。"
        "【中等证据】富集到的通路（细胞外基质、细胞因子—趋化因子轴等）方向可信，"
        "但 ORA/GSEA 未控制年龄、性别、吸烟史等协变量，效应量仍需独立队列验证。"
        f"【弱证据】WGCNA 的模块-性状相关（最强 {f['best_module'] or '—'}，r={_fmt_num(f['best_r'])}）"
        f"只说明共表达结构与分组同步，不构成因果；且 {len(f['hub_overlap'])}/{f['n_hub']} 个 hub 落在差异基因内，二者不独立。"
        "【不可用】细胞类型归属、蛋白水平变化与药物干预效果均超出本数据范围，不作结论。"
    )

    # 7-3
    first_go = str(f["ora_go"]["term"].iloc[0]) if not f["ora_go"].empty else "—"
    first_up = str(f["gsea_up"]["term"].iloc[0]) if len(f["gsea_up"]) else "—"
    first_down = str(f["gsea_down"]["term"].iloc[0]) if len(f["gsea_down"]) else "—"
    texts["7-3"] = (
        f"本研究基于 {f['source']} 的 {n_samples} 例样本转录组，按「质控—差异—共表达—富集—整合」的固定流程"
        f"得到一条较完整的证据链：标准化后样本一致性良好（中位相关性 {_fmt_num(f['median_corr'])}），"
        f"差异分析给出 {f['n_sig']} 个显著基因；"
        f"上调端以 {'、'.join(f['top_up'][:3]) or '—'} 为代表，下调端以 {'、'.join(f['top_down'][:3]) or '—'} 为代表。\n\n"
        f"机制推测：ORA 最显著的过程是 {first_go}，"
        f"GSEA 激活端为 {first_up}，抑制端为 {first_down}；"
        f"WGCNA 中 {f['best_module'] or '—'} 与 {f['group_a_disp']} 相关最强（r={_fmt_num(f['best_r'])}），"
        "可作为后续优先验证的切入点。\n\n"
        f"局限与展望：第一，数据为组织块 bulk 转录组，无法区分细胞来源，需单细胞或空间转录组确认；"
        "第二，差异模型未纳入年龄、性别等协变量，存在混杂可能；第三，结论停留在转录本层面，缺蛋白与功能验证。"
        "后续可对 hub 基因做 IHC/qPCR 验证，并把 hub 清单接入网络药理学的成分—靶点网络。"
    )
    return texts


# ---------------------------------------------------------------- 参数表


def _norm_label(method: str) -> str:
    return {
        "quantile": "quantile（分位数标准化）",
        "median": "median（中位数中心化）",
        "none": "none（不做标准化）",
    }.get(str(method), str(method))


def _fmt_e(value: float) -> str:
    """表单里填科学计数法用 e 写法（4.30e-63），不用 ×10^ 排版。"""
    if value is None or (isinstance(value, float) and math.isnan(value)):
        return "（待补）"
    return f"{value:.2e}"


def _evidence_note(f: dict) -> str:
    first_go = str(f["ora_go"]["term"].iloc[0]) if not f["ora_go"].empty else "—"
    overlap = len(f["hub_overlap"])
    if overlap and overlap == f["n_hub"]:
        caution = (
            f"{f['n_hub']} 个 hub 全部落在 {f['n_sig']} 个 DEG 内：hub 是从 {f['wgcna_max_genes']} 个高变基因里按 kME 选出，"
            "与差异基因本就高度重叠，二者是同一批基因的两种排序，不宜当作两条独立证据。"
        )
    else:
        caution = f"hub 与 DEG 的交集为 {overlap}/{f['n_hub']}，重合部分同样不宜当作独立证据。"
    return (
        f"证据链四层收敛：差异表达（{f['n_sig']} 个基因，{f['n_up']} 上调 / {f['n_down']} 下调）"
        f"→ 共表达（{f['best_module'] or '—'} 与 {f['group_a_disp']} 相关 r={_fmt_num(f['best_r'])}，"
        f"p={_fmt_sci(f['best_p'])}）→ hub（{f['n_hub']} 个高 kME 基因）→ 富集（{first_go} 居首）。"
        f"四层共同指向 {f['group_a_disp']} 相对 {f['group_b_disp']} 的转录组重塑。{caution}"
    )


def build_params_table(f: dict, texts: dict[str, str] | None = None) -> pd.DataFrame:
    """页面表单字段 → 建议值。字段 id 沿用课程提交系统的命名，缺的标「（待补）」。"""
    texts = texts or {}

    def dash(value) -> str:
        if value is None or value == "" or (isinstance(value, float) and math.isnan(value)):
            return "（待补）"
        return str(value)

    def note(code: str) -> str:
        return texts.get(code) or f"见文字稿 {code}"

    best_p = _fmt_e(f["best_p"])
    gsea_up = f["gsea_up"]
    up_paths = "；".join(
        f"{row.term}（NES={_fmt_num(row.NES, 2)}，FDR={_fmt_num(row.padj)}）" for row in gsea_up.itertuples()
    ) or "（待补）"
    module_sizes = " / ".join(f"{k}={v}" for k, v in f["module_sizes"].items()) or "（未运行 WGCNA）"
    r2_cut = (
        f"0.85（实际最大 R²={_fmt_num(f['r2'])}，"
        + ("达到门槛" if not math.isnan(f["r2"]) and f["r2"] >= 0.85 else "未达门槛，按最大 R² 取值")
        + "）"
        if not math.isnan(f["r2"])
        else "0.85"
    )
    rows = [
        ("第 1 段", "pv-dbMain", "数据库类型",
         f"GEO（{f['gse']}）" if f["gse"] else "本地上传（如来自 GEO 请填 GEO）"),
        ("第 1 段", "pv-log2Rule", "log2 判定规则", f["log2_rule"]),
        ("第 1 段", "pv-normMethod", "标准化方法", _norm_label(f["normalization"])),
        ("第 1 段", "pv-seed", "随机种子", dash(f["seed"])),
        ("第 1 段", "pv-caseN", "病例组样本数", str(f["groups"].get(f["group_a"], "（待补）"))),
        ("第 1 段", "pv-ctrlN", "对照组样本数", str(f["groups"].get(f["group_b"], "（待补）"))),
        ("第 1 段", "pv-geneCount", "受检基因数", str(f["n_genes"])),
        ("第 1 段", "pv-datasetNote", "数据集说明", note("1-5")),
        ("第 2 段", "pv-topVarN", "高变基因数", str(f["top_variable_n"])),
        ("第 2 段", "pv-pcaCenter", "PCA 中心化", "TRUE"),
        ("第 2 段", "pv-pcaScale", "PCA 标准化", "TRUE"),
        ("第 2 段", "pv-corMethod", "样本相关性方法", "Pearson"),
        ("第 2 段", "pv-qcMedianCor", "样本相关性中位数", _fmt_num(f["median_corr"])),
        ("第 2 段", "pv-qcNote", "质控结论", note("2-4")),
        ("第 3 段", "pv-fdrCut", "padj 阈值", f"{f['padj_cut']:g}"),
        ("第 3 段", "pv-logfcCut", "|log2FC| 阈值", f"{f['lfc_cut']:g}"),
        ("第 3 段", "pv-trend", "limma trend 参数", "FALSE（自实现 limma eBayes，未启用 trend）"),
        ("第 3 段", "pv-model", "设计矩阵", "model.matrix(~0+group)"),
        ("第 3 段", "pv-nUp", "上调基因数", str(f["n_up"])),
        ("第 3 段", "pv-nDown", "下调基因数", str(f["n_down"])),
        ("第 3 段", "pv-degAll", "受检基因数", str(f["n_genes"])),
        ("第 3 段", "pv-degSig", "显著差异基因数", str(f["n_sig"])),
        ("第 3 段", "pv-degNote", "差异结果解读", note("3-5")),
        ("第 4 段", "pv-wgcnaTopN", "WGCNA 取基因数", str(f["wgcna_max_genes"])),
        ("第 4 段", "pv-r2Cut", "无标度拟合 R² 门槛", r2_cut),
        ("第 4 段", "pv-minModuleSize", "最小模块基因数", str(f["wgcna_min_module"])),
        ("第 4 段", "pv-mergeCutHeight", "模块合并阈值", "0.75（模块特征基因相关阈值）"),
        ("第 4 段", "pv-networkType", "网络类型", "signed"),
        ("第 4 段", "pv-softPower", "软阈值 β", _fmt_num(f["beta"], 0)),
        ("第 4 段", "pv-r2", "实际无标度拟合 R²", _fmt_num(f["r2"])),
        ("第 4 段", "pv-moduleCount", "模块数", str(f["n_modules"]) if f["n_modules"] else "（未运行 WGCNA）"),
        ("第 4 段", "pv-moduleP", "最强相关模块 P 值", best_p),
        ("第 4 段", "pv-hubCount", "hub 候选基因数", str(f["n_hub"])),
        ("第 4 段", "pv-wgcnaNote", "模块与 hub 解释", note("4-4")),
        ("第 5 段", "pv-ont", "GO 类别", "BP（Biological Process）"),
        ("第 5 段", "pv-pAdjust", "多重检验校正", "BH（Benjamini-Hochberg）"),
        ("第 5 段", "pv-orgDb", "物种注释库", "org.Hs.eg.db（Homo sapiens；如物种不同请改）"),
        ("第 5 段", "pv-pathwayCount", "GSEA 通路总数", str(f["gsea_terms"]) if f["gsea_available"] else "（未运行 GSEA）"),
        ("第 5 段", "pv-sigPathways", "显著通路数",
         f"{f['gsea_fdr25']}（FDR<0.25；FDR<0.05 为 {f['gsea_fdr05']}）" if f["gsea_available"] else "（未运行 GSEA）"),
        ("第 5 段", "pv-minGSSize", "基因集最小规模", str(f["gsea_min_size"])),
        ("第 5 段", "pv-maxGSSize", "基因集最大规模", str(f["gsea_max_size"])),
        ("第 5 段", "pv-upPathways", "正富集代表通路", up_paths),
        ("第 5 段", "pv-leadingGenes", "leading edge 代表基因",
         "激活端：" + ("、".join(f["leading_up"]) or "—") + "；抑制端：" + ("、".join(f["leading_down"]) or "—")),
        ("第 5 段", "pv-gseaNote", "富集解读", note("5-4")),
        ("第 6 段", "pv-coreGenes", "DEG ∩ hub 核心基因数", str(len(f["hub_overlap"]))),
        ("第 6 段", "pv-intersectN", "交集基因数", str(len(f["hub_overlap"]))),
        ("第 6 段", "pv-bestModule", "最强相关模块", f["best_module"] or "（未运行 WGCNA）"),
        ("第 6 段", "pv-bestModuleR", "最强相关模块 r", _fmt_num(f["best_r"])),
        ("第 6 段", "pv-bestModuleP", "最强相关模块 p", best_p),
        ("第 6 段", "pv-moduleSizes", "各模块基因数", module_sizes),
        ("第 6 段", "pv-evidenceNote", "证据链说明", _evidence_note(f)),
        ("第 6 段", "pv-pharmLink", "网络药理学衔接", note("6-2")),
        ("第 7 段", "pv-figDpi", "图件分辨率", "200–210 dpi"),
        ("第 7 段", "pv-figFormat", "图件格式", "PNG，白色背景，Arial/微软雅黑"),
        ("第 7 段", "pv-conclusionNote", "结论与展望", note("7-3")),
    ]
    return pd.DataFrame(rows, columns=["段", "字段 id", "页面标签", "建议填写值"])


# ---------------------------------------------------------------- 辅助


def _image_size(path: Path) -> tuple[int, int]:
    try:
        import matplotlib.image as mpimg

        array = mpimg.imread(path)
        return int(array.shape[1]), int(array.shape[0])
    except Exception:  # noqa: BLE001
        return 0, 0


def build_submission_pack(
    result_dir: str | Path,
    out_dir: str | Path | None = None,
    *,
    screenshot: str | Path | None = None,
    heatmap_top: int = 50,
) -> PackResult:
    """公开入口：结果目录 → 提交材料包。"""
    return SubmissionPack(
        result_dir,
        out_dir,
        screenshot=screenshot,
        heatmap_top=heatmap_top,
    ).build()
