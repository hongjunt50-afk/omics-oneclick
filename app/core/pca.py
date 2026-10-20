"""主成分分析：标准化 + SVD，输出得分、载荷与解释方差。"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass
class PCAResult:
    scores: pd.DataFrame
    loadings: pd.DataFrame
    explained_variance: np.ndarray
    explained_ratio: np.ndarray
    n_components: int

    def to_frames(self) -> dict[str, pd.DataFrame]:
        ev = pd.DataFrame(
            {
                "component": [f"PC{i + 1}" for i in range(self.n_components)],
                "explained_variance": self.explained_variance,
                "explained_ratio": self.explained_ratio,
                "cumulative_ratio": np.cumsum(self.explained_ratio),
            }
        )
        return {"scores": self.scores, "loadings": self.loadings, "variance": ev}


def run_pca(
    expr: pd.DataFrame,
    *,
    n_components: int = 10,
    scale: bool = True,
    groups: pd.Series | None = None,
) -> PCAResult:
    """样本为观测、基因为变量。scale=True 时按基因做 z-score（推荐）。"""
    if expr.shape[1] < 2:
        raise ValueError("样本数不足 2 个，无法做 PCA")
    if expr.shape[0] < 2:
        raise ValueError("基因数不足，无法做 PCA")

    matrix = expr.to_numpy(dtype=float).T  # 样本 × 基因
    matrix = matrix - matrix.mean(axis=0, keepdims=True)
    if scale:
        sd = matrix.std(axis=0, ddof=1, keepdims=True)
        sd[sd == 0] = 1.0
        matrix = matrix / sd

    matrix = np.nan_to_num(matrix, nan=0.0, posinf=0.0, neginf=0.0)

    max_comp = min(matrix.shape[0], matrix.shape[1])
    k = int(min(max(n_components, 2), max_comp))

    u, s, vt = np.linalg.svd(matrix, full_matrices=False)
    variance = (s**2) / max(matrix.shape[0] - 1, 1)
    total = variance.sum()
    ratio = variance / total if total > 0 else np.zeros_like(variance)

    scores = u[:, :k] * s[:k]
    scores = pd.DataFrame(
        scores,
        index=[str(c) for c in expr.columns],
        columns=[f"PC{i + 1}" for i in range(k)],
    )
    loadings = pd.DataFrame(
        vt[:k].T,
        index=[str(g) for g in expr.index],
        columns=[f"PC{i + 1}" for i in range(k)],
    )

    if groups is not None:
        scores.insert(0, "group", groups.reindex(expr.columns).astype(str).values)

    return PCAResult(
        scores=scores,
        loadings=loadings,
        explained_variance=variance[:k],
        explained_ratio=ratio[:k],
        n_components=k,
    )
