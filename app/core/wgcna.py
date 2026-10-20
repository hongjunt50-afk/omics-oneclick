"""WGCNA：软阈值选择、邻接矩阵、TOM 相异度、层次聚类模块识别与模块-性状关联。

模块识别采用「TOM 层次聚类 + 分支切割 + 模块特征基因合并」的工程化实现，
与 R 包 WGCNA 的 dynamicTreeCut 不完全等价，但参数含义一致
（软阈值 β、模块最小基因数、模块合并阈值）。
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy import stats
from scipy.cluster.hierarchy import fcluster, leaves_list, linkage
from scipy.spatial.distance import squareform


@dataclass
class WGCNAResult:
    beta: float
    power_table: pd.DataFrame
    modules: pd.Series
    module_eigengenes: pd.DataFrame
    module_trait: pd.DataFrame
    traits: pd.DataFrame
    linkage: np.ndarray
    dissimilarity: np.ndarray
    gene_order: list[str]
    notes: list[str] = field(default_factory=list)

    def module_summary(self) -> pd.DataFrame:
        counts = self.modules.value_counts()
        rows = [{"module": m, "n_genes": int(n)} for m, n in counts.items() if m != "grey"]
        return pd.DataFrame(rows).sort_values("n_genes", ascending=False) if rows else pd.DataFrame(
            columns=["module", "n_genes"]
        )

    def hub_genes(self, module: str, expr: pd.DataFrame, top: int = 20) -> pd.DataFrame:
        """模块内按与模块特征基因的相关性排序（kME 的近似）。"""
        members = self.modules[self.modules == module].index
        if module not in self.module_eigengenes.columns or len(members) == 0:
            return pd.DataFrame(columns=["gene", "kME"])
        me = self.module_eigengenes[module]
        block = expr.loc[members]
        corr = np.array([_pearson_safe(block.loc[g].to_numpy(float), me.to_numpy(float)) for g in members])
        table = pd.DataFrame({"gene": list(members), "kME": corr})
        return table.reindex(table["kME"].abs().sort_values(ascending=False).index).head(top)


def _pearson_safe(a: np.ndarray, b: np.ndarray) -> float:
    if a.std() == 0 or b.std() == 0:
        return 0.0
    value = float(np.corrcoef(a, b)[0, 1])
    return 0.0 if not np.isfinite(value) else value


def _adjacency(corr: np.ndarray, beta: float, signed: bool) -> np.ndarray:
    a = ((1.0 + corr) / 2.0) ** beta if signed else np.abs(corr) ** beta
    np.fill_diagonal(a, 1.0)
    return a


def _corr_matrix(x: np.ndarray) -> np.ndarray:
    sd = x.std(axis=1, ddof=1, keepdims=True)
    sd[sd == 0] = 1.0
    z = (x - x.mean(axis=1, keepdims=True)) / sd
    corr = (z @ z.T) / max(x.shape[1] - 1, 1)
    corr = np.nan_to_num(corr)
    np.fill_diagonal(corr, 1.0)
    return np.clip(corr, -1.0, 1.0)


def scale_free_fit(adj: np.ndarray) -> float:
    """连接度分布 log-log 线性拟合；返回带符号 R²（标准做法符号取反，正值即合格）。"""
    k = adj.sum(axis=1) - 1.0
    k = k[k > 1e-6]
    if k.size < 10:
        return 0.0
    bins = int(min(20, max(8, k.size // 10)))
    counts, edges = np.histogram(k, bins=bins)
    centers = (edges[:-1] + edges[1:]) / 2
    valid = counts > 0
    if valid.sum() < 4:
        return 0.0
    lx = np.log10(centers[valid])
    ly = np.log10(counts[valid] / counts.sum())
    slope, intercept = np.polyfit(lx, ly, 1)
    pred = slope * lx + intercept
    ss_res = float(((ly - pred) ** 2).sum())
    ss_tot = float(((ly - ly.mean()) ** 2).sum())
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return float(-np.sign(slope) * max(r2, 0.0))


def pick_soft_threshold(
    corr: np.ndarray,
    *,
    powers: tuple[int, ...] = tuple(range(1, 21)),
    target_r2: float = 0.80,
    signed: bool = True,
) -> tuple[float, pd.DataFrame]:
    rows = []
    best_beta = float(powers[4] if len(powers) > 4 else powers[0])
    for beta in powers:
        adj = _adjacency(corr, beta, signed)
        r2 = scale_free_fit(adj)
        rows.append({"power": beta, "scale_free_R2": r2, "mean_connectivity": float((adj.sum(axis=1) - 1).mean())})
        if r2 >= target_r2:
            best_beta = float(beta)
            break
    else:
        finite = [r for r in rows if np.isfinite(r["scale_free_R2"])]
        if finite:
            best_beta = float(max(finite, key=lambda r: r["scale_free_R2"])["power"])

    table = pd.DataFrame(rows)
    table["selected"] = table["power"] == best_beta
    return best_beta, table


def _tom(adj: np.ndarray) -> np.ndarray:
    k = adj.sum(axis=1) - 1.0
    shared = adj @ adj
    denom = np.minimum.outer(k, k) + 1.0 - adj
    denom = np.where(denom <= 0, 1e-12, denom)
    tom = (shared + adj) / denom
    np.fill_diagonal(tom, 1.0)
    return np.clip(tom, 0.0, 1.0)


def _module_eigengenes(x: np.ndarray, genes: list[str], samples: list[str], modules: pd.Series) -> pd.DataFrame:
    out: dict[str, np.ndarray] = {}
    pos_of = {g: i for i, g in enumerate(genes)}
    for module in pd.unique(modules):
        if module == "grey":
            continue
        members = modules[modules == module].index
        pos = [pos_of[g] for g in members]
        block = x[pos].T
        block = block - block.mean(axis=0, keepdims=True)
        if block.shape[1] == 1:
            out[module] = block[:, 0]
            continue
        u, s, _ = np.linalg.svd(block, full_matrices=False)
        pc1 = u[:, 0] * s[0]
        if _pearson_safe(block.mean(axis=1), pc1) < 0:
            pc1 = -pc1
        out[module] = pc1
    return pd.DataFrame(out, index=samples)


def _merge_modules(modules: pd.Series, eigengenes: pd.DataFrame, *, min_size: int, merge_cut: float) -> pd.Series:
    modules = modules.copy()
    eig = eigengenes.copy()

    # 1) 迭代合并特征基因高度相关的模块
    while eig.shape[1] > 1:
        names = list(eig.columns)
        corr = np.corrcoef(eig.to_numpy().T)
        iu = np.triu_indices(len(names), k=1)
        values = corr[iu]
        best = int(np.argmax(values))
        if not np.isfinite(values[best]) or values[best] < merge_cut:
            break
        a, b = names[iu[0][best]], names[iu[1][best]]
        keep, drop = (a, b) if eig[a].var() >= eig[b].var() else (b, a)
        modules[modules == drop] = keep
        eig[keep] = eig[[keep, drop]].mean(axis=1)
        eig = eig.drop(columns=[drop])

    # 2) 小模块并入最相近的大模块，否则归 grey
    counts = modules.value_counts()
    big = [m for m in counts.index if m != "grey" and counts[m] >= min_size]
    for module, size in counts.items():
        if module == "grey" or size >= min_size or module not in eig.columns:
            continue
        scores = {
            big_mod: abs(_pearson_safe(eig[module].to_numpy(float), eig[big_mod].to_numpy(float)))
            for big_mod in big
            if big_mod in eig.columns
        }
        target = max(scores, key=scores.get) if scores else None
        if target is not None and scores[target] >= 0.3:
            modules[modules == module] = target
        else:
            modules[modules == module] = "grey"

    return modules


def run_wgcna(
    expr: pd.DataFrame,
    *,
    traits: pd.DataFrame | None = None,
    max_genes: int = 3000,
    signed: bool = True,
    beta: float | None = None,
    min_module_size: int = 30,
    merge_cut: float = 0.75,
    cut_height: float = 0.75,
) -> WGCNAResult:
    """expr 为 基因×样本；traits 为 样本×性状（数值型），用于模块-性状关联。"""
    if expr.shape[0] < 30:
        raise ValueError("基因数少于 30，无法进行 WGCNA（建议 ≥1000）")

    notes: list[str] = []
    variances = expr.var(axis=1).sort_values(ascending=False)
    kept = [g for g in variances.index if variances[g] > 0][:max_genes]
    if len(kept) < variances.size:
        notes.append(f"按方差取前 {len(kept)} 个基因（原 {variances.size} 个）用于 WGCNA")

    sub = expr.loc[kept]
    genes = [str(g) for g in sub.index]
    samples = [str(c) for c in sub.columns]
    sub.index = genes
    x = sub.to_numpy(dtype=float)

    corr = _corr_matrix(x)
    auto_beta, power_table = pick_soft_threshold(corr, signed=signed)
    chosen = float(beta) if beta is not None else auto_beta
    power_table["selected"] = power_table["power"] == chosen

    adj = _adjacency(corr, chosen, signed)
    tom = _tom(adj)
    diss = 1.0 - tom
    np.fill_diagonal(diss, 0.0)
    diss = (diss + diss.T) / 2.0

    tree = linkage(squareform(diss, checks=False), method="average")
    labels = fcluster(tree, t=cut_height, criterion="distance")
    modules = pd.Series([f"M{lab}" for lab in labels], index=genes)

    eig = _module_eigengenes(x, genes, samples, modules)
    modules = _merge_modules(modules, eig, min_size=min_module_size, merge_cut=merge_cut)
    eig = _module_eigengenes(x, genes, samples, modules)

    if traits is not None and not traits.empty:
        traits = traits.reindex(samples)
        module_trait = _module_trait_correlation(eig, traits)
    else:
        traits = pd.DataFrame()
        module_trait = pd.DataFrame()

    notes.append(f"软阈值 β = {chosen:g}，识别模块 {len([m for m in modules.unique() if m != 'grey'])} 个")

    order = [genes[i] for i in leaves_list(tree)]
    return WGCNAResult(
        beta=chosen,
        power_table=power_table,
        modules=modules,
        module_eigengenes=eig,
        module_trait=module_trait,
        traits=traits,
        linkage=tree,
        dissimilarity=diss,
        gene_order=order,
        notes=notes,
    )


def _module_trait_correlation(eigengenes: pd.DataFrame, traits: pd.DataFrame) -> pd.DataFrame:
    numeric = traits.apply(pd.to_numeric, errors="coerce")
    rows = []
    for module in eigengenes.columns:
        me = eigengenes[module].to_numpy(float)
        for trait in numeric.columns:
            y = numeric[trait].to_numpy(float)
            if np.isnan(y).all() or np.nanstd(y) == 0:
                continue
            r, p = stats.pearsonr(me, y)
            rows.append({"module": module, "trait": str(trait), "cor": float(r), "pvalue": float(p)})
    return pd.DataFrame(rows)
