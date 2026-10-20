"""数据载入、格式识别与分组解析。"""

from __future__ import annotations

import gzip
import io
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

SUPPORTED_SUFFIXES = {".csv", ".tsv", ".txt", ".tab", ".xlsx", ".xls", ".gz"}


@dataclass
class Dataset:
    """规整后的表达矩阵：expr 的行为基因，列为样本。"""

    expr: pd.DataFrame
    annotation: pd.DataFrame | None = None
    groups: pd.Series | None = None
    source: str = ""
    log2_transformed: bool = False
    notes: list[str] = field(default_factory=list)

    @property
    def n_genes(self) -> int:
        return int(self.expr.shape[0])

    @property
    def n_samples(self) -> int:
        return int(self.expr.shape[1])

    def group_names(self) -> list[str]:
        if self.groups is None:
            return []
        return [str(g) for g in pd.unique(self.groups)]

    def summary(self) -> str:
        parts = [f"{self.n_genes} 基因 × {self.n_samples} 样本"]
        if self.groups is not None:
            counts = self.groups.value_counts()
            detail = "、".join(f"{k}({v})" for k, v in counts.items())
            parts.append(f"分组：{detail}")
        return "；".join(parts)


def _strip_index(name: str) -> str:
    """去掉样本名尾部的编号与分隔符：chp_26 -> chp，ipf_H03 -> ipf_H。"""
    stem = re.sub(r"[\s_.\-]*\d+\s*$", "", name).strip(" _.-")
    return stem or name


def _infer_groups_detail(sample_names: list[str]) -> tuple[pd.Series | None, list[str]]:
    """按样本名推断分组，返回 (分组 Series, 说明)。无法推断时 Series 为 None。"""
    names = [str(n).strip() for n in sample_names]
    if not names:
        return None, []

    stems = {name: _strip_index(name) for name in names}
    counts: dict[str, int] = {}
    for stem in stems.values():
        counts[stem] = counts.get(stem, 0) + 1

    notes: list[str] = []
    # 单例词干若能在分隔符边界上截短到某个多例词干，就归并过去（ipf_H03 -> ipf）
    for name, stem in list(stems.items()):
        if counts.get(stem, 0) > 1:
            continue
        for match in reversed(list(re.finditer(r"[\s_.\-]+", stem))):
            prefix = stem[: match.start()]
            if prefix and counts.get(prefix, 0) > 1:
                stems[name] = prefix
                counts[stem] -= 1
                counts[prefix] += 1
                notes.append(f"样本 {name} 无标准编号，按前缀归入「{prefix}」")
                break

    series = pd.Series(stems)
    if series.nunique() < 2:
        return None, notes
    return series, notes


def infer_groups(sample_names: list[str]) -> pd.Series | None:
    """从样本名推断分组，例如 Control_1 / Treat_2 -> Control / Treat。

    单个样本名不合规则不再让整份推断作废，只有全部推断出同一组时才返回 None。
    """
    series, _notes = _infer_groups_detail(sample_names)
    return series


def group_table_frame(sample_names: list[str]) -> tuple[pd.DataFrame | None, list[str]]:
    """生成两列分组表（sample, group）及推断说明，供界面「自动生成分组表」使用。"""
    series, notes = _infer_groups_detail(sample_names)
    if series is None:
        return None, notes
    frame = pd.DataFrame({"sample": series.index.astype(str), "group": series.to_numpy()})
    frame = frame.sort_values(["group", "sample"], kind="stable").reset_index(drop=True)
    return frame, notes


def save_group_table(frame: pd.DataFrame, path: str | Path) -> Path:
    """写出分组表 CSV（UTF-8 BOM，Excel 双击不乱码）。"""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, encoding="utf-8-sig")
    return path


def _read_delimited(text: str, sep: str, name: str) -> pd.DataFrame:
    return pd.read_csv(io.StringIO(text), sep=sep, low_memory=False)


def _read_raw_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()
    opener = gzip.open if suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8-sig", errors="replace") as handle:
        text = handle.read()

    inner = path.name[:-3].lower() if suffix == ".gz" else path.name.lower()
    if inner.endswith((".xlsx", ".xls")):
        engine = "openpyxl" if inner.endswith(".xlsx") else "xlrd"
        with opener(path, "rb") as handle:
            return pd.read_excel(io.BytesIO(handle.read()), sheet_name=0, engine=engine)

    # 嗅探分隔符：优先制表符，其次逗号，再次分号
    first_line = next((line for line in text.splitlines() if line.strip()), "")
    for sep in ("\t", ",", ";"):
        if first_line.count(sep) >= 1:
            return _read_delimited(text, sep, path.name)
    return _read_delimited(text, r"\s+", path.name)


def _looks_numeric(series: pd.Series) -> bool:
    converted = pd.to_numeric(series, errors="coerce")
    return converted.notna().mean() >= 0.8


def _parse_annotation(df: pd.DataFrame, index_col: str) -> pd.DataFrame | None:
    """从主表里切出 *_meta 之类的元数据列；没有则返回 None。"""
    meta_cols = [c for c in df.columns if isinstance(c, str) and c.endswith("_meta")]
    if not meta_cols:
        return None
    ann = df.set_index(index_col)[meta_cols].copy()
    ann.columns = [c[: -len("_meta")] for c in ann.columns]
    return ann


def load_dataset(
    path: str | Path,
    *,
    genes_in_rows: bool | None = None,
    group_file: str | Path | None = None,
    force_log2: bool | None = None,
) -> Dataset:
    """读取表达矩阵并规整为 基因×样本。genes_in_rows=None 时自动判断方向。"""
    path = Path(path)
    if not path.exists():
        raise FileNotFoundError(f"文件不存在：{path}")

    raw = _read_raw_table(path)
    if raw.empty:
        raise ValueError("文件为空或无法解析出数据表")

    raw = raw.dropna(axis=1, how="all")
    if raw.shape[1] < 2:
        raise ValueError("数据表至少需要 1 个标识列 + 1 个样本列")

    notes: list[str] = []
    id_col = raw.columns[0]
    annotation = _parse_annotation(raw, id_col)

    data_cols = [c for c in raw.columns[1:] if c not in (annotation.columns if annotation is not None else [])]
    data_cols = [c for c in data_cols if not str(c).endswith("_meta")]
    body = raw[[id_col] + data_cols].copy()

    numeric_ratio = sum(_looks_numeric(body[c]) for c in data_cols) / max(len(data_cols), 1)
    if numeric_ratio < 0.5:
        # 首行其实是样本名（矩阵无表头）
        body = body.copy()
        new_cols = [str(id_col)] + [str(v) for v in body.iloc[0, 1:]]
        body.columns = new_cols
        body = body.iloc[1:].reset_index(drop=True)
        id_col = new_cols[0]
        notes.append("检测到无表头矩阵，已用首行作为样本名")

    body = body.set_index(id_col)

    if genes_in_rows is None:
        genes_in_rows = _detect_orientation(body)

    if genes_in_rows:
        expr = body.apply(pd.to_numeric, errors="coerce")
    else:
        expr = body.apply(pd.to_numeric, errors="coerce").T
        notes.append("检测到样本为行，已转置为 基因×样本")

    expr.index = [str(i) for i in expr.index]
    expr = expr[~expr.index.isna()]
    expr = expr.dropna(axis=0, how="all")
    expr = expr.fillna(0.0)

    duplicated = expr.index.duplicated()
    if duplicated.any():
        expr = expr.groupby(level=0, sort=False).mean()
        notes.append(f"合并了 {int(duplicated.sum())} 个重复基因名（取均值）")

    keep = expr.std(axis=1, ddof=0) > 0
    if (~keep).sum() > 0:
        notes.append(f"剔除 {int((~keep).sum())} 个零方差基因")
        expr = expr.loc[keep]

    if expr.shape[1] < 2:
        raise ValueError("有效样本数不足 2 个，无法分析")

    transformed = _apply_log(expr, force_log2)
    if transformed:
        expr = np.log2(expr + 1.0)
        notes.append("已做 log2(x+1) 变换")

    groups = None
    if group_file:
        groups = load_group_file(group_file, list(expr.columns))
        if groups is not None:
            notes.append(f"分组来自文件 {Path(group_file).name}")
    if groups is None:
        groups, infer_notes = _infer_groups_detail(list(expr.columns))
        if groups is not None:
            notes.append("按样本名推断分组")
            notes.extend(infer_notes)

    return Dataset(
        expr=expr,
        annotation=annotation,
        groups=groups,
        source=str(path),
        log2_transformed=bool(transformed),
        notes=notes,
    )


def _detect_orientation(body: pd.DataFrame) -> bool:
    """返回 True 表示基因在行。靠“行数远多于列数”和标识列是否像基因名判断。"""
    n_rows, n_cols = body.shape
    if n_rows >= n_cols * 2:
        return True
    if n_cols >= n_rows * 2:
        return False
    first = body.index.astype(str)
    sample_like = first.str.contains(r"\d+$", regex=True).mean()
    return sample_like < 0.5


def _apply_log(expr: pd.DataFrame, force_log2: bool | None) -> bool:
    if force_log2 is not None:
        return force_log2
    finite = expr.to_numpy(dtype=float)
    finite = finite[np.isfinite(finite)]
    if finite.size == 0:
        return False
    if (finite < 0).any():
        return False
    return float(np.percentile(finite, 99)) > 50.0


def load_group_file(path: str | Path, samples: list[str]) -> pd.Series | None:
    """分组表：两列（样本名、组名）或单列且以样本名为索引的表格。"""
    path = Path(path)
    df = _read_raw_table(path)
    if df.shape[1] == 1:
        series = df.iloc[:, 0]
        series.index = [str(c) for c in series.index]
    else:
        labels = df.columns[:2]
        series = pd.Series(df[labels[1]].astype(str).values, index=df[labels[0]].astype(str).values)

    series = series.reindex(samples)
    missing = series.isna().sum()
    if missing == len(samples):
        return None
    if missing:
        raise ValueError(f"分组表缺少 {missing} 个样本的组别")
    return series.astype(str)
