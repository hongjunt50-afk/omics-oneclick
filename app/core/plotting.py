"""统一视觉风格的绘图：火山图、聚类热图、PCA 散点、富集条形/气泡图、GSEA 曲线、WGCNA 树。"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import LinearSegmentedColormap  # noqa: E402

plt.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False
plt.rcParams["figure.dpi"] = 150
plt.rcParams["savefig.bbox"] = "tight"
plt.rcParams["savefig.facecolor"] = "white"

PALETTE = {
    "up": "#E64B5B",
    "down": "#3C8DBC",
    "ns": "#B9BFC7",
    "accent": "#2BB3A3",
    "accent2": "#7B61FF",
    "grid": "#E3E7EC",
    "text": "#20262E",
}

DIVERGING = LinearSegmentedColormap.from_list(
    "omics_div", ["#1F5FA9", "#5C9AD5", "#F5F6F8", "#EFA08A", "#C2304A"]
)
SEQUENTIAL = LinearSegmentedColormap.from_list("omics_seq", ["#F4F8FB", "#8FC8D8", "#2B7FA8", "#123D5C"])

GROUP_COLORS = ["#2BB3A3", "#E64B5B", "#5C6BC0", "#F2A03D", "#7B61FF", "#4DB6AC", "#D4587A", "#8D6E63"]


def _save(fig, path: Path, *, close: bool = True) -> Path:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=200, facecolor="white")
    if close:
        plt.close(fig)
    return path


def volcano_plot(
    deg_table: pd.DataFrame,
    path: Path,
    *,
    padj_cut: float = 0.05,
    lfc_cut: float = 1.0,
    title: str = "差异基因火山图",
    label_top: int = 12,
) -> Path:
    table = deg_table.copy()
    y = -np.log10(table["pvalue"].clip(lower=1e-300))
    x = table["log2FC"]

    sig_up = (table["padj"] < padj_cut) & (x >= lfc_cut)
    sig_down = (table["padj"] < padj_cut) & (x <= -lfc_cut)
    ns = ~(sig_up | sig_down)

    fig, ax = plt.subplots(figsize=(7.2, 6.0))
    ax.scatter(x[ns], y[ns], s=8, c=PALETTE["ns"], alpha=0.55, linewidths=0, label=f"不显著 ({int(ns.sum())})")
    ax.scatter(x[sig_up], y[sig_up], s=14, c=PALETTE["up"], alpha=0.85, linewidths=0, label=f"上调 ({int(sig_up.sum())})")
    ax.scatter(x[sig_down], y[sig_down], s=14, c=PALETTE["down"], alpha=0.85, linewidths=0, label=f"下调 ({int(sig_down.sum())})")

    ax.axhline(-np.log10(padj_cut), color="#8A94A6", ls="--", lw=0.9)
    ax.axvline(lfc_cut, color="#8A94A6", ls="--", lw=0.9)
    ax.axvline(-lfc_cut, color="#8A94A6", ls="--", lw=0.9)

    top = table[(sig_up | sig_down)].nsmallest(label_top, "pvalue")
    for _, row in top.iterrows():
        ax.annotate(
            row["gene"],
            (row["log2FC"], -np.log10(max(row["pvalue"], 1e-300))),
            fontsize=7,
            color=PALETTE["text"],
            xytext=(3, 3),
            textcoords="offset points",
        )

    ax.set_xlabel("log2 倍数变化 (log2FC)")
    ax.set_ylabel("-log10 P 值")
    ax.set_title(title, fontsize=13, color=PALETTE["text"], pad=12)
    ax.legend(loc="upper right", fontsize=8, frameon=False)
    ax.grid(True, color=PALETTE["grid"], lw=0.6, alpha=0.9)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    return _save(fig, path)


def heatmap_plot(
    expr: pd.DataFrame,
    genes: list[str],
    path: Path,
    *,
    groups: pd.Series | None = None,
    title: str = "差异基因聚类热图",
    max_genes: int = 60,
    cluster_samples: bool = True,
    stats: pd.DataFrame | None = None,
) -> Path:
    """差异基因表达热图。

    ``stats`` 传入含 gene / log2FC 的差异表时，输出「热图 + log2FC 条形图 + 组均值折线图」
    三联版式；不传则只画热图本体。
    """
    selected = [g for g in genes if g in expr.index][:max_genes]
    if not selected:
        raise ValueError("热图没有可用的差异基因")

    block = expr.loc[selected]
    z = block.sub(block.mean(axis=1), axis=0).div(block.std(axis=1).replace(0, 1), axis=0)
    z = z.clip(-3, 3)

    samples = list(z.columns)
    if cluster_samples and len(samples) > 2:
        order = _sample_order(z, samples)
        z = z[order]

    group_labels = None
    if groups is not None:
        group_labels = [str(groups.get(s, "")) for s in z.columns]
        unique = list(dict.fromkeys(group_labels))
        color_map = {g: GROUP_COLORS[i % len(GROUP_COLORS)] for i, g in enumerate(unique)}
    else:
        color_map = {}
    pretty = {"ipf": "IPF", "control": "Control", "ctrl": "Control", "chp": "CHP"}

    fc_map: dict[str, float] = {}
    if stats is not None and "gene" in stats.columns:
        fc_col = "log2FC" if "log2FC" in stats.columns else ("logFC" if "logFC" in stats.columns else None)
        if fc_col:
            for g, v in zip(stats["gene"].astype(str), stats[fc_col]):
                if pd.notna(v):
                    fc_map[g] = float(v)
    triple = bool(fc_map)
    if triple:  # 行按 log2FC 降序，热图呈上红下蓝，条形图与折线单调可读
        z = z.iloc[sorted(range(len(z.index)), key=lambda i: -fc_map.get(z.index[i], 0.0))]

    if triple:
        fig = plt.figure(figsize=(12.4, 9.6))
        grid = fig.add_gridspec(
            3,
            3,
            width_ratios=[1, 1, 0.026],
            height_ratios=[0.05, 1.66, 1.0],
            wspace=0.26,
            hspace=0.40,
            left=0.142,
            right=0.952,
            top=0.938,
            bottom=0.072,
        )
        group_ax = fig.add_subplot(grid[0, 0:2])
        ax = fig.add_subplot(grid[1, 0:2])
        cbar_ax = fig.add_subplot(grid[1:, 2])
        ax_b = fig.add_subplot(grid[2, 0])
        ax_c = fig.add_subplot(grid[2, 1])
    else:
        fig = plt.figure(
            figsize=(min(max(6.5, 0.34 * z.shape[1] + 3.6), 12.5), max(5.0, 0.16 * z.shape[0] + 3.0))
        )
        grid = fig.add_gridspec(
            2,
            2,
            width_ratios=[1, 0.028],
            height_ratios=[0.045, 1],
            wspace=0.02,
            hspace=0.05,
            left=0.145,
            right=0.955,
            top=0.925,
            bottom=0.06,
        )
        group_ax = fig.add_subplot(grid[0, 0])
        ax = fig.add_subplot(grid[1, 0])
        cbar_ax = fig.add_subplot(grid[1, 1])
        ax_b = ax_c = None

    im = ax.imshow(z.to_numpy(), aspect="auto", cmap=DIVERGING, vmin=-3, vmax=3, interpolation="nearest")
    ax.set_xticks(range(z.shape[1]))
    if z.shape[1] <= 40:
        ax.set_xticklabels(z.columns, rotation=90, fontsize=7)
    else:  # 样本多时只留分组色带，样本名省略
        ax.set_xticklabels([""] * z.shape[1])
    ax.set_yticks(range(z.shape[0]))
    ax.set_yticklabels(z.index, fontsize=6.5)
    # 标题画在画布层，避免与上方分组色带抢位置
    fig.text(0.02, 0.965, title, fontsize=12.5, ha="left", color=PALETTE["text"])
    ax.set_ylabel("基因", fontsize=9, labelpad=6)
    if z.shape[1] > 40:
        ax.set_xlabel("样本（按层次聚类排序，上方色带为分组）", fontsize=9, labelpad=6)
    for spine in ax.spines.values():
        spine.set_visible(False)
    if triple:
        ax.set_title("A　Top50 基因表达 z-score 热图（每行按 log2FC 降序）", fontsize=10.5, loc="left", pad=6)

    cbar = fig.colorbar(im, cax=cbar_ax)
    cbar.set_label("z-score", fontsize=8)
    cbar.ax.tick_params(labelsize=7)

    group_ax.set_xlim(-0.5, z.shape[1] - 0.5)
    group_ax.set_ylim(0, 1)
    group_ax.axis("off")
    if group_labels:
        for i, label in enumerate(group_labels):
            group_ax.add_patch(
                plt.Rectangle((i - 0.5, 0.1), 1, 0.8, color=color_map[label], lw=0)
            )
        handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in color_map.values()]
        group_ax.legend(
            handles,
            [pretty.get(str(k).lower(), str(k)) for k in color_map.keys()],
            loc="lower right",
            bbox_to_anchor=(1.0, 1.15),
            ncol=min(len(color_map), 8),
            fontsize=8,
            frameon=False,
            handlelength=1.1,
            handleheight=0.9,
            columnspacing=1.0,
            borderaxespad=0.0,
        )

    # ---------- B：各基因 log2FC 水平条形图（行序与热图一致） ----------
    if triple and ax_b is not None:
        fc = np.array([fc_map.get(g, np.nan) for g in z.index], dtype=float)
        ypos = np.arange(len(fc))
        colors = [PALETTE["up"] if v > 0 else PALETTE["down"] for v in np.nan_to_num(fc)]
        ax_b.barh(ypos, np.nan_to_num(fc), color=colors, height=0.7, edgecolor="none", alpha=0.9)
        ax_b.axvline(0, color="#6B7684", lw=0.9)
        ax_b.set_ylim(len(fc) - 0.5, -0.5)
        step = max(1, len(fc) // 10)
        show = list(range(0, len(fc), step))
        ax_b.set_yticks(show)
        ax_b.set_yticklabels([z.index[i] for i in show], fontsize=6)
        ax_b.set_xlabel("log2FC（IPF vs Control）", fontsize=8.5)
        ax_b.set_title("B　各基因的差异倍数", fontsize=10.5, loc="left", pad=6, color=PALETTE["text"])
        ax_b.grid(True, axis="x", color=PALETTE["grid"], lw=0.6)
        ax_b.set_axisbelow(True)
        for side in ("top", "right"):
            ax_b.spines[side].set_visible(False)

    # ---------- C：两组平均 z-score 折线 ----------
    if triple and ax_c is not None and group_labels:
        xs = np.arange(1, z.shape[0] + 1)
        key_order = list(color_map.keys())
        for key in key_order:
            mask = np.array([lab == key for lab in group_labels])
            if not mask.any():
                continue
            ax_c.plot(
                xs,
                z.iloc[:, mask].mean(axis=1).to_numpy(),
                marker="o",
                markersize=2.8,
                linewidth=1.4,
                color=color_map[key],
                label=pretty.get(str(key).lower(), str(key)),
            )
        ax_c.axhline(0, color="#6B7684", lw=0.8, ls="--", alpha=0.8)
        ax_c.set_xlabel("基因排名（按 |log2FC| 降序）", fontsize=8.5)
        ax_c.set_ylabel("组内平均 z-score", fontsize=8.5)
        ax_c.set_title("C　两组平均表达沿基因的变化", fontsize=10.5, loc="left", pad=6, color=PALETTE["text"])
        ax_c.legend(fontsize=8, frameon=False, loc="best")
        ax_c.grid(True, axis="y", color=PALETTE["grid"], lw=0.6)
        ax_c.set_axisbelow(True)
        for side in ("top", "right"):
            ax_c.spines[side].set_visible(False)

    return _save(fig, path)


def _sample_order(z: pd.DataFrame, samples: list[str]) -> list[str]:
    from scipy.cluster.hierarchy import leaves_list, linkage
    from scipy.spatial.distance import pdist

    if len(samples) < 3:
        return samples
    try:
        tree = linkage(pdist(z.to_numpy().T, metric="correlation"), method="average")
        return [samples[i] for i in leaves_list(tree)]
    except Exception:
        return samples


def pca_plot(
    scores: pd.DataFrame,
    explained: np.ndarray,
    path: Path,
    *,
    groups: pd.Series | None = None,
    title: str = "主成分分析 (PCA)",
) -> Path:
    fig, ax = plt.subplots(figsize=(6.6, 5.4))
    labels = groups.reindex(scores.index).astype(str) if groups is not None else pd.Series(
        ["样本"] * len(scores), index=scores.index
    )
    unique = list(dict.fromkeys(labels.tolist()))

    for i, group in enumerate(unique):
        mask = labels == group
        ax.scatter(
            scores.loc[mask, "PC1"],
            scores.loc[mask, "PC2"],
            s=70,
            alpha=0.85,
            color=GROUP_COLORS[i % len(GROUP_COLORS)],
            label=f"{group} (n={int(mask.sum())})",
            edgecolors="white",
            linewidths=0.8,
        )
        if mask.sum() > 2:
            ellipse = _confidence_ellipse(scores.loc[mask, "PC1"], scores.loc[mask, "PC2"])
            if ellipse is not None:
                ax.add_patch(ellipse)

    for name, row in scores.iterrows():
        ax.annotate(str(name), (row["PC1"], row["PC2"]), fontsize=6, alpha=0.7, xytext=(3, 3), textcoords="offset points")

    ax.axhline(0, color=PALETTE["grid"], lw=0.8)
    ax.axvline(0, color=PALETTE["grid"], lw=0.8)
    ax.set_xlabel(f"PC1 ({explained[0] * 100:.1f}%)")
    ax.set_ylabel(f"PC2 ({explained[1] * 100:.1f}%)")
    ax.set_title(title, fontsize=13, pad=12)
    ax.legend(fontsize=8, frameon=False)
    ax.grid(True, color=PALETTE["grid"], lw=0.6, alpha=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    return _save(fig, path)


def _confidence_ellipse(x: pd.Series, y: pd.Series):
    from matplotlib.patches import Ellipse

    if len(x) < 3:
        return None
    cov = np.cov(x, y)
    if not np.all(np.isfinite(cov)):
        return None
    values, vectors = np.linalg.eigh(cov)
    order = values.argsort()[::-1]
    values, vectors = values[order], vectors[:, order]
    if values[0] <= 0:
        return None
    angle = np.degrees(np.arctan2(vectors[1, 0], vectors[0, 0]))
    width, height = 2 * 2.0 * np.sqrt(values)
    return Ellipse(
        (float(x.mean()), float(y.mean())),
        width=width,
        height=height,
        facecolor="none",
        edgecolor="#9AA5B1",
        ls="--",
        lw=1.0,
        alpha=0.9,
    )


def enrichment_barplot(
    table: pd.DataFrame,
    path: Path,
    *,
    top: int = 15,
    title: str = "富集分析",
    direction: str | None = None,
) -> Path:
    data = table.copy()
    if direction == "up":
        data = data[data["direction"] == "activated"] if "direction" in data.columns else data
    elif direction == "down":
        data = data[data["direction"] == "suppressed"] if "direction" in data.columns else data

    value_col = "NES" if "NES" in data.columns else "fold_enrichment"
    data = data.sort_values(value_col, ascending=False).head(top)
    data = data.iloc[::-1]

    height = max(4.0, 0.36 * len(data) + 1.4)
    fig, ax = plt.subplots(figsize=(8.2, height))
    colors = [PALETTE["up"] if v >= 0 else PALETTE["down"] for v in data[value_col]]
    bars = ax.barh(range(len(data)), data[value_col], color=colors, alpha=0.85, height=0.72)

    ax.set_yticks(range(len(data)))
    labels = [_shorten(t, 62) for t in data["term"]]
    ax.set_yticklabels(labels, fontsize=7.5)
    ax.set_xlabel("NES（归一化富集分数）" if value_col == "NES" else "富集倍数 (fold enrichment)")
    ax.set_title(title, fontsize=12, pad=10)
    ax.axvline(0, color="#8A94A6", lw=0.9)

    for bar, value in zip(bars, data[value_col]):
        ax.text(
            value + (0.02 * np.nanmax(np.abs(data[value_col])) if value >= 0 else -0.02 * np.nanmax(np.abs(data[value_col]))),
            bar.get_y() + bar.get_height() / 2,
            f"{value:.2f}",
            va="center",
            ha="left" if value >= 0 else "right",
            fontsize=6.5,
            color=PALETTE["text"],
        )

    ax.grid(True, axis="x", color=PALETTE["grid"], lw=0.6, alpha=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    return _save(fig, path)


def enrichment_bubble(
    table: pd.DataFrame,
    path: Path,
    *,
    top: int = 20,
    title: str = "富集气泡图",
) -> Path:
    data = table.copy()
    pcol = "padj" if "padj" in data.columns else "pvalue"
    data = data.nsmallest(top, pcol).iloc[::-1]
    xcol = "fold_enrichment" if "fold_enrichment" in data.columns else "NES"
    sizecol = "overlap" if "overlap" in data.columns else "set_size"

    fig, ax = plt.subplots(figsize=(8.2, max(4.0, 0.34 * len(data) + 1.6)))
    score = -np.log10(data[pcol].clip(lower=1e-300))
    sizes = 40 + 320 * (data[sizecol] - data[sizecol].min()) / max(data[sizecol].max() - data[sizecol].min(), 1)
    scatter = ax.scatter(
        data[xcol], range(len(data)), s=sizes, c=score, cmap=SEQUENTIAL, alpha=0.9, edgecolors="white", linewidths=0.8
    )
    ax.set_yticks(range(len(data)))
    ax.set_yticklabels([_shorten(t, 60) for t in data["term"]], fontsize=7.5)
    ax.set_xlabel("富集倍数" if xcol == "fold_enrichment" else "NES")
    ax.set_title(title, fontsize=12, pad=10)
    ax.grid(True, axis="x", color=PALETTE["grid"], lw=0.6, alpha=0.8)
    ax.set_axisbelow(True)
    cbar = fig.colorbar(scatter, ax=ax, pad=0.02, fraction=0.035)
    cbar.set_label("-log10 P", fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    return _save(fig, path)


def gsea_plot(
    metric: pd.Series,
    hits: list[str],
    path: Path,
    *,
    term: str,
    nes: float,
    pvalue: float,
    title: str | None = None,
) -> Path:
    genes = list(metric.index)
    values = metric.to_numpy(dtype=float)
    n = len(genes)
    where = {g: i for i, g in enumerate(genes)}
    positions = np.array(sorted(where[g] for g in hits if g in where), dtype=np.int64)
    if positions.size == 0:
        raise ValueError("基因集与排序列表无交集")

    work = np.abs(values) + 1e-9
    w = work[positions]
    cum = np.cumsum(w) / w.sum()
    miss = (positions - np.arange(positions.size)) / max(n - positions.size, 1)
    running = cum - miss
    xs = np.arange(n)

    fig, (ax_top, ax_bottom) = plt.subplots(
        2, 1, figsize=(7.6, 5.2), sharex=True, gridspec_kw={"height_ratios": [2.0, 1.0], "hspace": 0.08}
    )

    curve = np.zeros(n)
    cursor = 0
    for i in range(n):
        while cursor < positions.size and positions[cursor] <= i:
            cursor += 1
        curve[i] = running[cursor - 1] if cursor > 0 else 0.0

    ax_top.plot(xs, curve, color=PALETTE["accent"] if nes >= 0 else PALETTE["down"], lw=1.6)
    ax_top.axhline(0, color="#8A94A6", lw=0.8)
    peak = int(np.argmax(np.abs(curve)))
    ax_top.axvline(peak, color="#C2304A", ls=":", lw=0.9)
    ax_top.set_ylabel("富集分数 ES")
    ax_top.set_title(
        title or f"{_shorten(term, 56)}\nNES = {nes:.2f}    P = {pvalue:.3g}", fontsize=10.5, pad=8
    )
    ax_top.grid(True, color=PALETTE["grid"], lw=0.6, alpha=0.75)
    ax_top.set_axisbelow(True)
    for side in ("top", "right"):
        ax_top.spines[side].set_visible(False)

    ax_bottom.vlines(positions, 0, 1, color="#3B4350", lw=0.6)
    ax_bottom.set_ylim(0, 1)
    ax_bottom.set_yticks([])
    ax_bottom.set_xlabel("基因排序（按排序指标由高到低）")
    for side in ("top", "right", "left"):
        ax_bottom.spines[side].set_visible(False)

    ax_top.set_xlim(0, n)
    return _save(fig, path)


def wgcna_dendrogram(result, path: Path, *, title: str = "WGCNA 基因聚类树与模块划分") -> Path:
    import sys

    from scipy.cluster.hierarchy import dendrogram

    fig, (ax_tree, ax_mod) = plt.subplots(
        2, 1, figsize=(10.0, 4.6), sharex=True, gridspec_kw={"height_ratios": [3.0, 0.8], "hspace": 0.05}
    )
    # scipy 的 dendrogram 是递归实现，3000 个基因的树会超过默认 1000 层限制
    limit = sys.getrecursionlimit()
    sys.setrecursionlimit(max(limit, 4 * len(result.linkage) + 1000))
    try:
        dendrogram(result.linkage, ax=ax_tree, no_labels=True, color_threshold=0.75, above_threshold_color="#7F8A99")
    finally:
        sys.setrecursionlimit(limit)
    ax_tree.set_ylabel("TOM 相异度")
    ax_tree.set_title(title, fontsize=12, pad=10)
    for side in ("top", "right", "bottom"):
        ax_tree.spines[side].set_visible(False)
    ax_tree.set_xticks([])

    modules = result.modules.reindex(result.gene_order)
    uniq = list(dict.fromkeys(modules.tolist()))
    color_of = {m: (GROUP_COLORS[i % len(GROUP_COLORS)] if m != "grey" else "#C9CED6") for i, m in enumerate(uniq)}
    for i, module in enumerate(modules):
        ax_mod.add_patch(plt.Rectangle((i * 10, 0), 10, 1, color=color_of[module], lw=0))
    ax_mod.set_xlim(0, len(modules) * 10)
    ax_mod.set_ylim(0, 1)
    ax_mod.set_yticks([])
    ax_mod.set_xlabel("基因（按聚类顺序）")
    for side in ("top", "right", "left"):
        ax_mod.spines[side].set_visible(False)
    return _save(fig, path)


def wgcna_power_plot(power_table: pd.DataFrame, path: Path) -> Path:
    fig, ax = plt.subplots(figsize=(6.2, 4.2))
    ax.plot(power_table["power"], power_table["scale_free_R2"], "o-", color=PALETTE["accent"], lw=1.6, ms=5)
    ax.axhline(0.8, color="#C2304A", ls="--", lw=1.0, label="目标 R² = 0.8")
    selected = power_table[power_table["selected"]]
    if not selected.empty:
        ax.scatter(
            selected["power"], selected["scale_free_R2"], s=90, facecolors="none", edgecolors="#C2304A", lw=1.6, zorder=5,
            label=f"选用 β = {selected['power'].iloc[0]:g}",
        )
    ax.set_xlabel("软阈值 β")
    ax.set_ylabel("无标度拓扑拟合 R²")
    ax.set_title("软阈值选择", fontsize=12, pad=10)
    ax.grid(True, color=PALETTE["grid"], lw=0.6, alpha=0.8)
    ax.set_axisbelow(True)
    ax.legend(fontsize=8, frameon=False)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    return _save(fig, path)


def module_trait_plot(module_trait: pd.DataFrame, path: Path) -> Path:
    pivot = module_trait.pivot(index="module", columns="trait", values="cor")
    fig, ax = plt.subplots(figsize=(max(4.2, 1.1 * pivot.shape[1] + 2.4), max(3.2, 0.42 * pivot.shape[0] + 1.8)))
    im = ax.imshow(pivot.to_numpy(), cmap=DIVERGING, vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(pivot.shape[1]))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=8)
    ax.set_yticks(range(pivot.shape[0]))
    ax.set_yticklabels(pivot.index, fontsize=8)
    for i in range(pivot.shape[0]):
        for j in range(pivot.shape[1]):
            value = pivot.to_numpy()[i, j]
            if np.isfinite(value):
                ax.text(j, i, f"{value:.2f}", ha="center", va="center", fontsize=7,
                        color="white" if abs(value) > 0.55 else PALETTE["text"])
    ax.set_title("模块与性状相关性", fontsize=12, pad=10)
    cbar = fig.colorbar(im, ax=ax, fraction=0.035, pad=0.02)
    cbar.set_label("Pearson r", fontsize=8)
    cbar.ax.tick_params(labelsize=7)
    return _save(fig, path)


def qc_boxplot(expr: pd.DataFrame, path: Path) -> Path:
    n = expr.shape[1]
    # 样本多时宽度封顶，否则 200+ 样本会生成上万像素的横幅图
    fig, ax = plt.subplots(figsize=(min(max(6.0, 0.42 * n + 2.8), 15.0), 4.3))
    data = [expr[c].to_numpy() for c in expr.columns]
    box = ax.boxplot(data, patch_artist=True, showfliers=False, widths=0.62)
    for i, patch in enumerate(box["boxes"]):
        patch.set_facecolor(GROUP_COLORS[i % len(GROUP_COLORS)])
        patch.set_alpha(0.55)
        patch.set_edgecolor("#4A5568")
    for median in box["medians"]:
        median.set_color("#20262E")
    ax.set_xticks(range(1, n + 1))
    if n <= 40:
        ax.set_xticklabels(expr.columns, rotation=90, fontsize=7)
    else:  # 密集时每 10 个样本标一次名，避免糊成一团
        ax.set_xticklabels(
            [str(c) if i % 10 == 0 else "" for i, c in enumerate(expr.columns)],
            rotation=90,
            fontsize=6,
        )
    ax.set_ylabel("表达值（log2 尺度优先）")
    ax.set_title("样本表达分布", fontsize=12, pad=10)
    ax.grid(True, axis="y", color=PALETTE["grid"], lw=0.6, alpha=0.8)
    ax.set_axisbelow(True)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    return _save(fig, path)


def _as_group_series(groups, columns) -> "pd.Series | None":
    """把多种形态的分组信息统一成 index=sample、value=group 的字符串 Series。"""
    if groups is None:
        return None
    if isinstance(groups, pd.DataFrame):
        frame = groups.copy()
        lookup = {str(c).lower(): c for c in frame.columns}
        if "group" in lookup:
            key = lookup.get("sample", frame.columns[0])
            frame = frame.set_index(key)
            frame.index = frame.index.astype(str)
            ser = frame[lookup["group"]]
        elif frame.shape[1] == 1:
            ser = frame.iloc[:, 0]
        else:
            return None
    elif isinstance(groups, pd.Series):
        ser = groups
    else:
        return None
    ser = ser.copy()
    ser.index = ser.index.astype(str)
    ser = ser.dropna().astype(str).reindex([str(c) for c in columns])
    return ser


def qc_correlation_plot(corr: pd.DataFrame, path: Path, groups=None) -> Path:
    """样本相关性：左主图矩阵热图 + 下排组内/组间箱线图与分布直方图。"""
    corr = corr.copy()
    corr.index = corr.index.astype(str)
    corr.columns = corr.columns.astype(str)
    mat = corr.to_numpy(dtype=float)
    n = corr.shape[0]
    series = _as_group_series(groups, corr.columns)
    uniq: list[str] = []
    if series is not None and int(series.notna().sum()) == n:
        uniq = list(dict.fromkeys(series.tolist()))

    fig = plt.figure(figsize=(10.4, 7.8))
    grid = fig.add_gridspec(
        2,
        2,
        width_ratios=[1.55, 1.0],
        height_ratios=[1.0, 1.0],
        left=0.062,
        right=0.975,
        top=0.868,
        bottom=0.098,
        wspace=0.24,
        hspace=0.42,
    )
    ax_h = fig.add_subplot(grid[:, 0])
    lo = float(np.nanmin(mat))
    lo = max(0.0, np.floor(lo * 20) / 20)
    im = ax_h.imshow(mat, cmap=SEQUENTIAL, vmin=lo, vmax=1.0, interpolation="nearest", aspect="auto")
    step = max(1, int(round(n / 8)))
    ticks = list(range(0, n, step))
    ax_h.set_xticks(ticks)
    ax_h.set_yticks(ticks)
    ax_h.set_xticklabels([str(t + 1) for t in ticks], fontsize=7)
    ax_h.set_yticklabels([str(t + 1) for t in ticks], fontsize=7)
    ax_h.set_xlabel("样本序号", fontsize=9, labelpad=4)
    ax_h.set_ylabel("样本序号", fontsize=9, labelpad=4)
    ax_h.set_title("A　样本间 Pearson 相关系数矩阵", fontsize=11.5, loc="left", pad=10, color=PALETTE["text"])
    for spine in ax_h.spines.values():
        spine.set_edgecolor("#9AA5B1")
        spine.set_linewidth(0.8)
    cbar = fig.colorbar(im, ax=ax_h, fraction=0.046, pad=0.03)
    cbar.set_label("Pearson r", fontsize=8.5)
    cbar.ax.tick_params(labelsize=7)
    cbar.outline.set_edgecolor("#9AA5B1")

    ax_b = fig.add_subplot(grid[0, 1])
    ax_c = fig.add_subplot(grid[1, 1])

    # ---------- B：组内 / 组间 相关系数分布 ----------
    display = {"ipf": "IPF", "control": "Control", "ctrl": "Control", "chp": "CHP"}
    palette_b = ["#2BB3A3", "#5C6BC0", "#F2A03D", "#7B61FF", "#E64B5B"]
    buckets: list[tuple[str, np.ndarray, str]] = []
    if len(uniq) >= 2:
        k = 0
        for i, ga in enumerate(uniq):
            for gb in uniq[i:]:
                ia = series.index[series == ga].tolist()
                ib = series.index[series == gb].tolist()
                if not ia or not ib:
                    continue
                da = display.get(ga.lower(), ga)
                db = display.get(gb.lower(), gb)
                if ga == gb:
                    sub = corr.loc[ia, ia].to_numpy(dtype=float)
                    vals = sub[np.triu_indices(len(ia), k=1)]
                    label = f"{da} 组内"
                    note = f"(n={len(ia)})"
                else:
                    vals = corr.loc[ia, ib].to_numpy(dtype=float).ravel()
                    label = f"{da} — {db}"
                    note = f"(n={len(ia) * len(ib)})"
                buckets.append((f"{label}\n{note}", np.asarray(vals, dtype=float), palette_b[k % len(palette_b)]))
                k += 1
    if not buckets:
        mask = ~np.eye(n, dtype=bool)
        buckets.append((f"全部样本对\n(n={int(mask.sum() / 2)})", mat[mask], "#2BB3A3"))

    bp = ax_b.boxplot(
        [b[1] for b in buckets],
        patch_artist=True,
        widths=0.55,
        showfliers=False,
        medianprops=dict(color=PALETTE["text"], lw=1.4),
        whiskerprops=dict(color="#6B7684", lw=0.9),
        capprops=dict(color="#6B7684", lw=0.9),
    )
    for patch, bucket in zip(bp["boxes"], buckets):
        patch.set_facecolor(bucket[2])
        patch.set_alpha(0.5)
        patch.set_edgecolor("#4A5568")
        patch.set_linewidth(0.8)
    ax_b.set_xticks(range(1, len(buckets) + 1))
    ax_b.set_xticklabels([b[0] for b in buckets], fontsize=7.5)
    ax_b.set_ylabel("Pearson r", fontsize=9)
    ax_b.set_title("B　组内与组间相关性分布", fontsize=11.5, loc="left", pad=8, color=PALETTE["text"])
    ax_b.grid(True, axis="y", color=PALETTE["grid"], lw=0.7)
    ax_b.set_axisbelow(True)
    for side in ("top", "right"):
        ax_b.spines[side].set_visible(False)

    # ---------- C：相关系数直方图 ----------
    for label, vals, color in buckets:
        ax_c.hist(
            vals,
            bins=60,
            histtype="stepfilled",
            alpha=0.45,
            color=color,
            edgecolor=color,
            linewidth=0.9,
            label=label.replace("\n", " "),
        )
    ax_c.set_xlabel("Pearson r", fontsize=9)
    ax_c.set_ylabel("样本对数量", fontsize=9)
    ax_c.set_title("C　相关系数频数分布", fontsize=11.5, loc="left", pad=8, color=PALETTE["text"])
    ax_c.legend(fontsize=7, frameon=False, loc="upper left")
    ax_c.grid(True, axis="y", color=PALETTE["grid"], lw=0.7)
    ax_c.set_axisbelow(True)
    for side in ("top", "right"):
        ax_c.spines[side].set_visible(False)

    fig.suptitle(f"样本相关性评估（n = {n} 例）", fontsize=13.5, x=0.062, ha="left", y=0.962, color=PALETTE["text"])
    fig.text(
        0.062,
        0.925,
        "对角线为样本自身（r = 1）；色块越深表示两个样本的整体表达谱越接近。组内 r 应高于组间 r，否则提示分组建模或批次效应问题。",
        fontsize=8.2,
        color="#5A6472",
        ha="left",
    )
    return _save(fig, path)


def _shorten(text: str, limit: int) -> str:
    text = str(text)
    return text if len(text) <= limit else text[: limit - 1] + "…"
