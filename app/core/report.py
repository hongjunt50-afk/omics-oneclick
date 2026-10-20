"""生成自包含 HTML 报告：图片内嵌 base64，表格可折叠，无外部依赖。"""

from __future__ import annotations

import base64
import html
import json
from datetime import datetime
from pathlib import Path

import pandas as pd

from .data_loader import Dataset
from .pipeline import AnalysisConfig, RunResult

ACCENT = "#2BB3A3"
ACCENT2 = "#3FA9F5"
UP = "#E64B5B"
DOWN = "#3C8DBC"

FIGURE_GROUPS = [
    ("质控与标准化", ("qc_",)),
    ("差异基因分析", ("volcano_", "heatmap_")),
    ("主成分分析", ("pca",)),
    ("富集分析 (ORA)", ("ora_",)),
    ("GSEA 分析", ("gsea",)),
    ("WGCNA 分析", ("wgcna",)),
]


def _img_tag(path: Path, caption: str) -> str:
    data = base64.b64encode(path.read_bytes()).decode("ascii")
    return (
        f'<figure><img src="data:image/png;base64,{data}" alt="{html.escape(caption)}"/>'
        f"<figcaption>{html.escape(caption)}</figcaption></figure>"
    )


def _table_html(frame: pd.DataFrame, max_rows: int = 30, float_fmt: str = "{:.4g}") -> str:
    if frame.empty:
        return '<p class="muted">（无数据）</p>'
    head = frame.head(max_rows)
    cols = list(head.columns)
    thead = "".join(f"<th>{html.escape(str(c))}</th>" for c in cols)
    body_rows = []
    for index, row in head.iterrows():
        cells = []
        if frame.index.name or index is not None:
            cells.append(f"<th>{html.escape(str(index))}</th>")
        for value in row:
            if isinstance(value, float):
                text = float_fmt.format(value)
            else:
                text = str(value)
            cells.append(f"<td>{html.escape(text)}</td>")
        body_rows.append("<tr>" + "".join(cells) + "</tr>")
    index_head = "<th></th>" if frame.index.name is not None else ""
    more = ""
    if len(frame) > max_rows:
        more = f'<p class="muted">仅显示前 {max_rows} 行，完整结果见 tables/ 目录下的 CSV。</p>'
    return (
        f'<div class="table-wrap"><table><thead><tr>{index_head}{thead}</tr></thead>'
        f"<tbody>{''.join(body_rows)}</tbody></table></div>{more}"
    )


def _read_csv(path: Path, rows: int | None = None) -> pd.DataFrame:
    kwargs = {"index_col": 0}
    if rows:
        kwargs["nrows"] = rows
    try:
        return pd.read_csv(path, encoding="utf-8-sig", **kwargs)
    except Exception:
        return pd.DataFrame()


def build_report(result: RunResult, dataset: Dataset, config: AnalysisConfig) -> Path:
    project_dir = result.project_dir
    summary = result.summary
    stamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

    # ---- 概览卡片 ----
    cards = [
        ("基因数", f"{dataset.n_genes:,}"),
        ("样本数", f"{dataset.n_samples:,}"),
        ("分析用时", f"{summary.get('elapsed_seconds', 0)} 秒"),
        ("产出图表", f"{len(result.figures)} 张"),
    ]
    if isinstance(summary.get("deg"), dict) and summary["deg"]:
        total_up = sum(v["n_up"] for v in summary["deg"].values())
        total_down = sum(v["n_down"] for v in summary["deg"].values())
        cards.append(("差异基因", f"↑{total_up} / ↓{total_down}"))
    if isinstance(summary.get("wgcna"), dict):
        cards.append(("WGCNA 模块", str(summary["wgcna"].get("modules", "-"))))

    card_html = "".join(
        f'<div class="card"><div class="card-label">{html.escape(k)}</div><div class="card-value">{html.escape(v)}</div></div>'
        for k, v in cards
    )

    # ---- 分组信息 ----
    group_html = '<p class="muted">未识别到分组信息（仅完成无监督分析）。</p>'
    if dataset.groups is not None:
        counts = dataset.groups.value_counts()
        rows = "".join(
            f"<tr><td>{html.escape(str(g))}</td><td>{int(n)}</td></tr>" for g, n in counts.items()
        )
        group_html = (
            '<div class="table-wrap narrow"><table><thead><tr><th>组别</th><th>样本数</th></tr></thead>'
            f"<tbody>{rows}</tbody></table></div>"
        )

    # ---- 数据说明 ----
    notes_html = "".join(f"<li>{html.escape(n)}</li>" for n in dataset.notes) or "<li>无</li>"

    # ---- 图组 ----
    sections = []
    used: set[str] = set()
    for title, prefixes in FIGURE_GROUPS:
        items = [
            (name, path)
            for name, path in result.figures.items()
            if name not in used and any(name.startswith(p) for p in prefixes)
        ]
        if not items:
            continue
        used.update(name for name, _ in items)
        figs = "".join(_img_tag(path, path.stem.replace("_", " ")) for _, path in items)
        sections.append(f'<section class="panel"><h2>{html.escape(title)}</h2><div class="grid">{figs}</div></section>')

    leftovers = [(n, p) for n, p in result.figures.items() if n not in used]
    if leftovers:
        figs = "".join(_img_tag(path, path.stem.replace("_", " ")) for _, path in leftovers)
        sections.append(f'<section class="panel"><h2>其他图表</h2><div class="grid">{figs}</div></section>')

    # ---- 关键表格 ----
    table_sections = []
    for name, path in result.tables.items():
        if name == "run_log" or not path.exists():
            continue
        if not name.startswith(("10_DEG", "30_ORA", "40_GSEA", "53_", "01_", "20_PCA_variance")):
            continue
        frame = _read_csv(path, rows=30)
        if frame.empty:
            continue
        label = (
            name.replace("10_DEG_", "差异分析 · ")
            .replace("30_ORA_", "富集分析 · ")
            .replace("40_GSEA_", "GSEA · ")
            .replace("53_WGCNA_module_summary", "WGCNA 模块汇总")
            .replace("01_sample_qc", "样本质控")
            .replace("20_PCA_variance", "PCA 方差解释")
        )
        table_sections.append(
            f'<details class="panel"><summary><h2>{html.escape(label)}</h2></summary>{_table_html(frame)}</details>'
        )

    # ---- 日志 ----
    log_html = "".join(
        f'<div class="log-line {("err" if "[ERROR]" in line else "warn" if "[WARNING]" in line else "")}">'
        f"{html.escape(line)}</div>"
        for line in result.log_lines
    )

    errors_html = ""
    if result.errors:
        items = "".join(f"<li>{html.escape(e)}</li>" for e in result.errors)
        errors_html = f'<section class="panel danger"><h2>未完成的步骤</h2><ul>{items}</ul></section>'

    doc = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="utf-8"/>
<meta name="viewport" content="width=device-width, initial-scale=1"/>
<title>OmicsOne 分析报告 · {html.escape(project_dir.name)}</title>
<style>
  :root {{
    --accent: {ACCENT}; --accent2: {ACCENT2}; --up: {UP}; --down: {DOWN};
    --ink: #17212B; --muted: #6B7785; --line: #E4E9EF; --bg: #F6F8FA;
  }}
  * {{ box-sizing: border-box; }}
  body {{ margin: 0; background: var(--bg); color: var(--ink);
         font: 15px/1.65 "Microsoft YaHei", "PingFang SC", system-ui, sans-serif; }}
  header {{ background: linear-gradient(120deg, #0E2030 0%, #123A47 55%, #0E2A33 100%); color: #EAF6F4;
            padding: 34px 44px 30px; }}
  header .brand {{ display: flex; align-items: center; gap: 12px; font-size: 13px; letter-spacing: .12em;
                   text-transform: uppercase; color: #8FD9CE; }}
  header h1 {{ margin: 10px 0 6px; font-size: 30px; font-weight: 600; letter-spacing: .01em; }}
  header .meta {{ color: #9FB6BD; font-size: 13px; }}
  main {{ max-width: 1180px; margin: 0 auto; padding: 26px 24px 60px; }}
  .cards {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(150px, 1fr)); gap: 14px; margin-bottom: 26px; }}
  .card {{ background: #fff; border: 1px solid var(--line); border-radius: 12px; padding: 16px 18px;
           border-left: 3px solid var(--accent); }}
  .card-label {{ font-size: 12px; color: var(--muted); letter-spacing: .04em; }}
  .card-value {{ font-size: 24px; font-weight: 600; margin-top: 6px; }}
  .panel {{ background: #fff; border: 1px solid var(--line); border-radius: 14px; padding: 20px 24px;
            margin-bottom: 22px; }}
  .panel h2 {{ font-size: 17px; margin: 0 0 14px; padding-left: 11px; border-left: 3px solid var(--accent); }}
  details.panel summary {{ cursor: pointer; list-style: none; }}
  details.panel summary::-webkit-details-marker {{ display: none; }}
  details.panel summary h2::after {{ content: "  ▾"; color: var(--accent); font-size: 13px; }}
  .grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(330px, 1fr)); gap: 18px; }}
  figure {{ margin: 0; background: #fff; border: 1px solid var(--line); border-radius: 10px; padding: 10px; }}
  figure img {{ width: 100%; height: auto; display: block; border-radius: 6px; }}
  figcaption {{ font-size: 12px; color: var(--muted); text-align: center; padding-top: 8px; }}
  .table-wrap {{ overflow-x: auto; }}
  .table-wrap.narrow {{ max-width: 420px; }}
  table {{ border-collapse: collapse; width: 100%; font-size: 13px; }}
  th, td {{ padding: 7px 10px; border-bottom: 1px solid var(--line); text-align: left; white-space: nowrap; }}
  thead th {{ background: #F0F5F7; font-weight: 600; color: #33505C; }}
  tbody tr:hover {{ background: #FAFCFD; }}
  .muted {{ color: var(--muted); font-size: 13px; }}
  ul {{ margin: 8px 0 0 18px; padding: 0; }}
  .log {{ max-height: 320px; overflow-y: auto; font: 12px/1.7 Consolas, "Courier New", monospace;
          background: #0F1B24; color: #CFE3E9; border-radius: 10px; padding: 14px 16px; }}
  .log-line {{ white-space: pre-wrap; }}
  .log-line.warn {{ color: #F4C77B; }}
  .log-line.err {{ color: #F08A8A; }}
  .panel.danger {{ border-color: #F3C6C6; background: #FFF8F8; }}
  footer {{ text-align: center; color: var(--muted); font-size: 12px; padding: 10px 0 0; }}
  @media print {{ body {{ background: #fff; }} .panel, .card {{ break-inside: avoid; }} }}
</style>
</head>
<body>
<header>
  <div class="brand">OmicsOne · 一键组学分析平台</div>
  <h1>组学分析报告</h1>
  <div class="meta">数据集：{html.escape(Path(dataset.source).name)}　　生成时间：{stamp}　　输出目录：{html.escape(str(project_dir))}</div>
</header>
<main>
  <div class="cards">{card_html}</div>

  <section class="panel">
    <h2>数据概况</h2>
    <p>{html.escape(dataset.summary())}</p>
    <h3 style="font-size:14px;margin:16px 0 8px;color:#33505C;">分组构成</h3>
    {group_html}
    <h3 style="font-size:14px;margin:16px 0 8px;color:#33505C;">数据处理说明</h3>
    <ul>{notes_html}</ul>
    <h3 style="font-size:14px;margin:16px 0 8px;color:#33505C;">分析参数</h3>
    <ul>
      <li>显著性阈值：padj &lt; {config.padj_cut}，|log2FC| ≥ {config.lfc_cut}</li>
      <li>标准化方法：{html.escape(config.normalization)}</li>
      <li>富集库：{html.escape("、".join(config.ora_libraries) or "未启用")}</li>
      <li>GSEA 库：{html.escape(config.gsea_library)}（置换 {config.gsea_permutations} 次）</li>
      <li>WGCNA：{"启用" if config.run_wgcna else "未启用"}{f"，基因上限 {config.wgcna_max_genes}" if config.run_wgcna else ""}</li>
    </ul>
  </section>

  {errors_html}
  {''.join(sections)}
  <section class="panel"><h2>结果表格</h2><p class="muted">点击展开查看；完整表格见 tables/ 目录。</p>
    {''.join(table_sections) or '<p class="muted">（无）</p>'}
  </section>
  <section class="panel"><h2>运行日志</h2><div class="log">{log_html}</div></section>
  <footer>OmicsOne v1.0 · 全部计算在本机完成，数据未离开本地</footer>
</main>
</body>
</html>"""

    report_path = project_dir / "分析报告.html"
    report_path.write_text(doc, encoding="utf-8")
    result.tables["report"] = report_path

    (project_dir / "分析概览.md").write_text(_markdown_summary(result, dataset, config), encoding="utf-8")
    return report_path


def _markdown_summary(result: RunResult, dataset: Dataset, config: AnalysisConfig) -> str:
    lines = [
        f"# OmicsOne 分析概览",
        "",
        f"- 数据集：`{Path(dataset.source).name}`",
        f"- 规模：{dataset.summary()}",
        f"- 输出目录：`{result.project_dir}`",
        f"- 用时：{result.summary.get('elapsed_seconds', 0)} 秒",
        "",
        "## 参数",
        "",
        f"- padj < {config.padj_cut}，|log2FC| ≥ {config.lfc_cut}",
        f"- 标准化：{config.normalization}",
        f"- 富集库：{', '.join(config.ora_libraries)}",
        f"- GSEA 库：{config.gsea_library}（{config.gsea_permutations} 次置换）",
        "",
        "## 产出",
        "",
        f"- 图表 {len([k for k in result.figures])} 张（figures/）",
        f"- 表格 {len([k for k in result.tables if k not in ('report', 'run_log')])} 个（tables/）",
        "",
        "## 目录结构",
        "",
        "```",
        f"{result.project_dir.name}/",
        "├── 分析报告.html        自包含报告（图片已内嵌）",
        "├── 分析概览.md          本文件",
        "├── summary.json        机器可读结果摘要",
        "├── run_log.txt         运行日志",
        "├── figures/            全部图表 PNG",
        "└── tables/             全部结果 CSV（UTF-8 BOM，Excel 可直接打开）",
        "```",
    ]
    if result.errors:
        lines += ["", "## 未完成的步骤", ""] + [f"- {e}" for e in result.errors]
    return "\n".join(lines) + "\n"
