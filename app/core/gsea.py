"""预排序 GSEA（preranked）：基因集置换法求 ES/NES/p/q 与核心基因。"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from .genesets import GeneSetLibrary


@dataclass
class GSEAResult:
    table: pd.DataFrame
    library: str
    n_genes: int
    n_permutations: int
    ranking: str

    def significant(self, *, padj: float = 0.25) -> pd.DataFrame:
        return self.table[self.table["padj"] < padj].copy()

    def leading_edge(self, term: str) -> list[str]:
        row = self.table[self.table["term"] == term]
        if row.empty:
            return []
        return [g for g in str(row.iloc[0]["leading_edge"]).split(",") if g]


def rank_metric(deg_table: pd.DataFrame, how: str = "log2FC") -> pd.Series:
    """按指定指标给基因排序，返回 基因 -> 排序值（已去重、降序）。"""
    genes = deg_table["gene"].astype(str)
    if how == "pvalue":
        p = pd.Series(-np.log10(deg_table["pvalue"].to_numpy(dtype=float).clip(min=1e-300)), index=genes)
        series = p * np.sign(deg_table["log2FC"].to_numpy(dtype=float))
    else:
        col = how if how in deg_table.columns else "log2FC"
        if col not in deg_table.columns:
            raise ValueError(f"缺少排序指标列：{how}")
        series = pd.Series(deg_table[col].to_numpy(dtype=float), index=genes)

    series = series[~series.index.duplicated()]
    series = series.replace([np.inf, -np.inf], np.nan).dropna()
    return series.sort_values(ascending=False)


def _running_sum(positions: np.ndarray, n_genes: int, w: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """返回命中位点处的跃迁后值 after、跃迁前值 before（形状同 positions）。"""
    h = positions.size
    n_miss = n_genes - h
    total_w = w.sum()
    cumw = np.cumsum(w)
    miss_frac = (positions - np.arange(h)) / n_miss
    after = cumw / total_w - miss_frac
    before = (cumw - w) / total_w - miss_frac
    return after, before


def _es_from_hits(positions: np.ndarray, n_genes: int, abs_rank: np.ndarray) -> tuple[float, np.ndarray, np.ndarray, np.ndarray]:
    w = abs_rank[positions]
    after, before = _running_sum(positions, n_genes, w)
    i_max = int(np.argmax(after))
    i_min = int(np.argmin(before))
    if after[i_max] >= -before[i_min]:
        return float(after[i_max]), positions[: i_max + 1], after, before
    return float(before[i_min]), positions[i_min:], after, before


def _null_es_by_size(
    size: int,
    n_genes: int,
    abs_rank: np.ndarray,
    n_perm: int,
    rng: np.random.Generator,
    block: int = 200,
) -> np.ndarray:
    """对固定基因集大小生成 n_perm 个随机集的 ES（全向量化）。"""
    out = np.empty(n_perm)
    filled = 0
    n_miss = n_genes - size
    if n_miss <= 0:
        return np.zeros(n_perm)
    while filled < n_perm:
        take = min(block, n_perm - filled)
        keys = rng.random((take, n_genes))
        picks = np.argpartition(keys, size - 1, axis=1)[:, :size]
        picks.sort(axis=1)

        w = abs_rank[picks]
        total = w.sum(axis=1)
        total[total <= 0] = 1.0
        cumw = np.cumsum(w, axis=1)
        miss_frac = (picks - np.arange(size)) / n_miss

        after = cumw / total[:, None] - miss_frac
        before = (cumw - w) / total[:, None] - miss_frac

        es_max = after.max(axis=1)
        es_min = before.min(axis=1)
        out[filled : filled + take] = np.where(es_max >= -es_min, es_max, es_min)
        filled += take
    return out


def preranked_gsea(
    metric: pd.Series,
    library: GeneSetLibrary,
    *,
    n_permutations: int = 1000,
    min_set_size: int = 15,
    max_set_size: int = 500,
    ranking: str = "log2FC",
    seed: int = 20260101,
    progress=None,
) -> GSEAResult:
    metric = metric.dropna().sort_values(ascending=False)
    genes = metric.index.to_numpy().astype(str)
    values = metric.to_numpy(dtype=float)
    n_genes = genes.size
    if n_genes < 100:
        raise ValueError("可排序基因数少于 100，GSEA 结果不可靠")

    abs_rank = np.maximum(np.abs(values), 1e-9)
    where = {g.upper(): i for i, g in enumerate(genes)}
    rng = np.random.default_rng(seed)

    candidate: list[tuple[str, np.ndarray]] = []
    for set_name, members in library.sets.items():
        upper = {g.upper() for g in members}
        if len(upper) > max_set_size:
            continue
        positions = np.array(sorted(where[g] for g in upper if g in where), dtype=np.int64)
        if positions.size < min_set_size or positions.size >= n_genes:
            continue
        candidate.append((set_name, positions))

    if not candidate:
        raise ValueError("没有符合大小要求的基因集，请换库或调整集合大小上下限")

    null_cache: dict[int, np.ndarray] = {}
    rows = []
    total_terms = len(candidate)

    for index, (set_name, positions) in enumerate(candidate):
        es, lead, after, before = _es_from_hits(positions, n_genes, abs_rank)
        h = positions.size

        null = null_cache.get(h)
        if null is None:
            null = _null_es_by_size(h, n_genes, abs_rank, n_permutations, rng)
            null_cache[h] = null

        if es >= 0:
            same = null[null >= 0]
            norm = float(same.mean()) if same.size else 1.0
            p = float((same >= es).mean()) if same.size else 1.0
        else:
            same = np.abs(null[null < 0])
            norm = float(same.mean()) if same.size else 1.0
            p = float((same >= abs(es)).mean()) if same.size else 1.0

        nes = es / norm if norm > 0 else 0.0

        rows.append(
            {
                "term": set_name,
                "set_size": h,
                "ES": es,
                "NES": nes,
                "pvalue": max(p, 1.0 / (n_permutations + 1)),
                "direction": "activated" if es >= 0 else "suppressed",
                "leading_edge": ",".join(genes[lead][:300]),
                "leading_edge_size": int(lead.size),
            }
        )
        if progress and (index % 200 == 0 or index == total_terms - 1):
            progress(index + 1, total_terms)

    table = pd.DataFrame(rows)
    table["padj"] = _gsea_fdr(table)
    table = table.sort_values("NES", ascending=False).reset_index(drop=True)
    return GSEAResult(
        table=table, library=library.name, n_genes=n_genes, n_permutations=n_permutations, ranking=ranking
    )


def _gsea_fdr(table: pd.DataFrame) -> np.ndarray:
    """对置换 p 值做 BH 校正（fgsea 的 padj 即此口径）。"""
    p = table["pvalue"].to_numpy(dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0.0, 1.0)
    return out
