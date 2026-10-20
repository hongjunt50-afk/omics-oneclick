"""差异基因分析：limma 风格的线性模型 + 经验贝叶斯调节 t 检验。"""

from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import numpy as np
import pandas as pd
from scipy import special, stats


@dataclass
class DEGResult:
    table: pd.DataFrame
    group_a: str
    group_b: str
    n_up: int
    n_down: int

    def significant(self, *, padj: float = 0.05, lfc: float = 1.0) -> pd.DataFrame:
        return self.table[(self.table["padj"] < padj) & (self.table["log2FC"].abs() >= lfc)].copy()


def _bh_adjust(pvalues: np.ndarray) -> np.ndarray:
    p = np.asarray(pvalues, dtype=float)
    n = p.size
    if n == 0:
        return p
    order = np.argsort(p)
    ranked = p[order] * n / np.arange(1, n + 1)
    ranked = np.minimum.accumulate(ranked[::-1])[::-1]
    out = np.empty(n, dtype=float)
    out[order] = np.clip(ranked, 0.0, 1.0)
    return out


def _fit_ebayes(residual_var: np.ndarray, df_residual: float) -> tuple[np.ndarray, float]:
    """矩估计法拟合逆卡方先验（limma 的 squeezeVar 精简版）。"""
    n = residual_var.size
    if n < 3 or df_residual <= 0:
        return residual_var, np.inf

    z = np.log(residual_var + 1e-12)
    e = z - special.digamma(df_residual / 2) + np.log(df_residual / 2)

    e_mean = float(np.mean(e))
    var_e = float(np.var(e, ddof=1))
    trigamma = float(special.polygamma(1, df_residual / 2))

    target = var_e - trigamma
    if target <= 0 or not np.isfinite(target):
        df_prior = np.inf
        var_prior = float(np.exp(e_mean))
    else:
        df_prior = max(_solve_df_prior(target), 0.1)
        var_prior = float(np.exp(e_mean + special.digamma(df_prior / 2) - np.log(df_prior / 2)))

    if not np.isfinite(var_prior) or var_prior <= 0:
        return residual_var, np.inf
    return var_prior, df_prior


def _solve_df_prior(target: float) -> float:
    """求 d0 使 trigamma(d0/2) = target，牛顿迭代。"""
    d0 = 4.0
    for _ in range(60):
        f = float(special.polygamma(1, d0 / 2)) - target
        fp = -float(special.polygamma(2, d0 / 2)) / 2
        if fp == 0:
            break
        d0_next = d0 - f / fp
        if d0_next <= 0.05 or not np.isfinite(d0_next):
            d0_next = d0 / 2
        if abs(d0_next - d0) < 1e-6:
            d0 = d0_next
            break
        d0 = d0_next
    return min(d0, 200.0)


def _squeeze(residual_var: np.ndarray, df_residual: float, var_prior: float, df_prior: float):
    if not np.isfinite(df_prior) or df_prior <= 0:
        return residual_var, np.full_like(residual_var, df_residual)
    post_var = (df_prior * var_prior + df_residual * residual_var) / (df_prior + df_residual)
    return post_var, np.full_like(residual_var, df_residual + df_prior)


def differential_expression(
    expr: pd.DataFrame,
    groups: pd.Series,
    group_a: str,
    group_b: str,
    *,
    paired: pd.Series | None = None,
) -> DEGResult:
    """group_a 相对 group_b 的差异。log2FC = mean(a) - mean(b)。"""
    samples = list(expr.columns)
    labels = groups.reindex(samples).astype(str)
    mask = labels.isin([group_a, group_b])
    if mask.sum() < 4:
        raise ValueError("两组样本数太少（至少各 2 个）")

    sub = expr.loc[:, mask]
    labels = labels[mask]

    n_a = int((labels == group_a).sum())
    n_b = int((labels == group_b).sum())
    if n_a < 2 or n_b < 2:
        raise ValueError(f"每组至少需要 2 个样本（当前 {group_a}={n_a}, {group_b}={n_b}）")

    values = sub.to_numpy(dtype=float)
    is_a = (labels == group_a).to_numpy()

    design = np.zeros((values.shape[1], 2))
    design[is_a, 0] = 1.0
    design[~is_a, 1] = 1.0

    if paired is not None:
        block = paired.reindex(sub.columns).astype(str).to_numpy()
        if len(set(block)) < len(block):
            dummies = pd.get_dummies(pd.Series(block), drop_first=True).to_numpy(dtype=float)
            # 配对块与分组若完全共线则报错
            stacked = np.hstack([design, dummies])
            if np.linalg.matrix_rank(stacked) < stacked.shape[1]:
                raise ValueError("配对信息与分组完全共线，无法估计")
            design = stacked

    coef, *_ = np.linalg.lstsq(design, values.T, rcond=None)
    fitted = (design @ coef).T
    residuals = values - fitted

    df_residual = values.shape[1] - design.shape[1]
    if df_residual <= 0:
        raise ValueError("样本数不足以估计模型（每组样本过少）")

    residual_var = (residuals**2).sum(axis=1) / df_residual
    residual_var = np.maximum(residual_var, 1e-12)

    var_prior, df_prior = _fit_ebayes(residual_var, df_residual)
    post_var, post_df = _squeeze(residual_var, df_residual, var_prior, df_prior)

    # 对照向量：group_a 系数 - group_b 系数
    contrast = np.zeros(design.shape[1])
    contrast[0] = 1.0
    contrast[1] = -1.0

    xtx_inv = np.linalg.pinv(design.T @ design)
    scale = float(np.sqrt(contrast @ xtx_inv @ contrast))

    estimate = coef.T @ contrast
    se = np.sqrt(post_var) * scale
    t_stat = estimate / np.maximum(se, 1e-12)

    p_values = 2 * stats.t.sf(np.abs(t_stat), post_df)
    padj = _bh_adjust(p_values)

    table = pd.DataFrame(
        {
            "gene": sub.index.astype(str),
            "log2FC": estimate,
            "mean_A": values[:, is_a].mean(axis=1),
            "mean_B": values[:, ~is_a].mean(axis=1),
            "t": t_stat,
            "pvalue": p_values,
            "padj": padj,
            "AveExpr": values.mean(axis=1),
        }
    )
    table["direction"] = np.where(table["log2FC"] >= 0, "up", "down")
    table = table.sort_values("pvalue").reset_index(drop=True)

    n_up = int(((table["padj"] < 0.05) & (table["log2FC"] >= 1)).sum())
    n_down = int(((table["padj"] < 0.05) & (table["log2FC"] <= -1)).sum())
    return DEGResult(table=table, group_a=str(group_a), group_b=str(group_b), n_up=n_up, n_down=n_down)


def pairwise_groups(groups: pd.Series) -> list[tuple[str, str]]:
    names = [str(g) for g in pd.unique(groups.astype(str))]
    return [(a, b) for a, b in combinations(names, 2)]
