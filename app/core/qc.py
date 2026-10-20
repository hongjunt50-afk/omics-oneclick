"""质控与标准化：样本统计、相关性、离群检测、log2 与分位数归一。"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd


@dataclass
class QCSummary:
    stats: pd.DataFrame
    correlation: pd.DataFrame
    outliers: list[str]
    notes: list[str]

    def to_frame(self) -> pd.DataFrame:
        return self.stats


def sample_qc(expr: pd.DataFrame, *, cor_threshold: float = 0.8) -> QCSummary:
    values = expr.to_numpy(dtype=float)
    stats = pd.DataFrame(
        {
            "sample": expr.columns.astype(str),
            "median": np.median(values, axis=0),
            "mean": values.mean(axis=0),
            "sd": values.std(axis=0, ddof=1),
            "iqr": np.percentile(values, 75, axis=0) - np.percentile(values, 25, axis=0),
            "min": values.min(axis=0),
            "max": values.max(axis=0),
            "detected": (values > 0).sum(axis=0),
        }
    )

    correlation = pd.DataFrame(
        np.corrcoef(values.T), index=expr.columns.astype(str), columns=expr.columns.astype(str)
    )

    notes: list[str] = []
    outliers: list[str] = []
    if correlation.shape[0] > 2:
        mask = ~np.eye(correlation.shape[0], dtype=bool)
        median_corr = np.where(mask, correlation.to_numpy(), np.nan)
        per_sample = np.nanmedian(median_corr, axis=1)
        overall = float(np.nanmedian(per_sample))
        for name, value in zip(correlation.index, per_sample):
            if value < overall - 0.15:
                outliers.append(str(name))
        if outliers:
            notes.append(f"{len(outliers)} 个样本与其他样本相关性明显偏低，建议复核：{'、'.join(outliers)}")
        else:
            notes.append(f"样本间中位相关性 {overall:.3f}，未发现明显离群样本")

    high = correlation.to_numpy()[~np.eye(correlation.shape[0], dtype=bool)]
    if high.size and high.max() >= 0.9999:
        notes.append("存在几乎完全相同的样本，请确认是否重复上样")

    return QCSummary(stats=stats, correlation=correlation, outliers=outliers, notes=notes)


def quantile_normalize(expr: pd.DataFrame) -> pd.DataFrame:
    """分位数标准化：让每个样本的分布对齐到平均经验分布。"""
    values = expr.to_numpy(dtype=float)
    ranks = np.argsort(np.argsort(values, axis=0), axis=0)
    sorted_values = np.sort(values, axis=0)
    reference = sorted_values.mean(axis=1)
    normalized = reference[ranks]
    return pd.DataFrame(normalized, index=expr.index, columns=expr.columns)


def median_center(expr: pd.DataFrame) -> pd.DataFrame:
    values = expr.to_numpy(dtype=float)
    centered = values - np.median(values, axis=0, keepdims=True)
    return pd.DataFrame(centered, index=expr.index, columns=expr.columns)


def standardize(expr: pd.DataFrame, method: str = "quantile") -> tuple[pd.DataFrame, str]:
    """返回 (标准化后矩阵, 说明文本)。method: none / quantile / median。"""
    if method in ("none", "", None):
        return expr, "未做样品间标准化"
    if method == "median":
        return median_center(expr), "已做中位数中心化"
    if method == "quantile":
        return quantile_normalize(expr), "已做分位数标准化（quantile normalization）"
    raise ValueError(f"未知标准化方法：{method}")
