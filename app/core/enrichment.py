"""过表达分析（ORA）：超几何检验 + BH 校正，背景取实测基因集。"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from .genesets import GeneSetLibrary


@dataclass
class ORAResult:
    table: pd.DataFrame
    library: str
    n_query: int
    n_background: int

    def significant(self, *, padj: float = 0.05) -> pd.DataFrame:
        return self.table[self.table["padj"] < padj].copy()


def over_representation(
    genes: list[str],
    background: list[str],
    library: GeneSetLibrary,
    *,
    min_set_size: int = 5,
    max_set_size: int = 500,
) -> ORAResult:
    query = {str(g).upper() for g in genes if str(g).strip()}
    universe = {str(g).upper() for g in background if str(g).strip()}

    if not query:
        raise ValueError("待富集的基因列表为空")
    if len(universe) < 10:
        raise ValueError("背景基因数过少（<10），无法做富集分析")

    query &= universe
    if len(query) < 3:
        raise ValueError("待富集基因与表达谱交集少于 3 个，无法做富集分析")

    n_query = len(query)
    n_bg = len(universe)
    rows = []

    for set_name, members in library.sets.items():
        member_set = {g.upper() for g in members} & universe
        size = len(member_set)
        if size < min_set_size or size > max_set_size:
            continue
        hits = member_set & query
        k = len(hits)
        if k < 2:
            continue
        p = float(stats.hypergeom.sf(k - 1, n_bg, size, n_query))
        expected = n_query * size / n_bg
        rows.append(
            {
                "term": set_name,
                "overlap": k,
                "set_size": size,
                "query_size": n_query,
                "expected": expected,
                "fold_enrichment": k / expected if expected > 0 else np.nan,
                "pvalue": p,
                "genes": ",".join(sorted(hits)[:80]),
            }
        )

    if not rows:
        raise ValueError(
            "没有任何基因集与待富集基因重叠 ≥2 个，请更换基因集库或放宽阈值"
        )

    table = pd.DataFrame(rows)
    table["padj"] = _bh(table["pvalue"].to_numpy())
    table = table.sort_values(["padj", "pvalue"]).reset_index(drop=True)
    return ORAResult(table=table, library=library.name, n_query=n_query, n_background=n_bg)


def _bh(p: np.ndarray) -> np.ndarray:
    p = np.asarray(p, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(ranked, 0, 1)
    return out


def deg_genes_for_ora(deg_table: pd.DataFrame, *, padj: float, lfc: float, direction: str = "both") -> list[str]:
    mask = deg_table["padj"] < padj
    if direction == "up":
        mask &= deg_table["log2FC"] > 0
    elif direction == "down":
        mask &= deg_table["log2FC"] < 0
    mask &= deg_table["log2FC"].abs() >= lfc
    return deg_table.loc[mask, "gene"].astype(str).tolist()
