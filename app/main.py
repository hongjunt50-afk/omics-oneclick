"""应用入口：设置高 DPI、字体与暗色主题后启动主窗口。"""

from __future__ import annotations

import sys
import traceback
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from PySide6.QtGui import QFont, QIcon
from PySide6.QtWidgets import QApplication, QMessageBox

from app.gui.main_window import MainWindow
from app.gui.theme import APP_QSS, logo_path


def _excepthook(exc_type, exc_value, exc_tb) -> None:
    """兜底：任何未捕获异常都弹窗而不是静默崩掉。"""
    text = "".join(traceback.format_exception(exc_type, exc_value, exc_tb))
    sys.stderr.write(text)
    try:
        QMessageBox.critical(None, "程序异常", f"发生未处理的错误：\n\n{exc_value}\n\n详细信息已写入控制台。")
    except Exception:
        pass


def _selftest(report_file: str) -> int:
    """无界面自检：验证打包后的资源、依赖与核心分析链路，结果写入 JSON。"""
    import json

    report: dict[str, object] = {"frozen": bool(getattr(sys, "_MEIPASS", None))}
    try:
        from app.core import genesets

        directory = genesets.asset_dir()
        report["asset_dir"] = str(directory)
        report["gmt_files"] = sorted(p.name for p in directory.glob("*.gmt"))
        libraries = genesets.available_libraries()
        report["libraries"] = {name: label for name, (_path, label) in libraries.items()}

        import numpy as np
        import pandas as pd

        from app.core.data_loader import Dataset
        from app.core.gsea import preranked_gsea
        from app.core.pipeline import AnalysisConfig, Runner

        rng = np.random.default_rng(3)
        kegg = genesets.load_library("KEGG_2021_Human")
        genes = sorted({g for member in kegg.sets.values() for g in member})[:1200]
        metric = pd.Series(rng.normal(0, 1, len(genes)), index=genes)
        metric.iloc[:40] += 3.0
        gsea = preranked_gsea(metric.sort_values(ascending=False), kegg, n_permutations=100)
        report["gsea_terms"] = len(gsea.table)

        samples = [f"Control_{i}" for i in range(1, 4)] + [f"Treat_{i}" for i in range(1, 4)]
        expr = pd.DataFrame(rng.normal(6, 1, (300, 6)), index=genes[:300], columns=samples)
        expr.iloc[:20, 3:] += 3.0
        dataset = Dataset(
            expr=expr,
            groups=pd.Series(["Control"] * 3 + ["Treat"] * 3, index=samples),
            source="selftest",
            log2_transformed=True,
        )
        target = Path(report_file).resolve().parent / "_selftest_run"
        config = AnalysisConfig(
            output_dir=target,
            group_a="Control",
            group_b="Treat",
            gsea_permutations=100,
            wgcna_max_genes=200,
            wgcna_min_module=10,
            ora_libraries=["KEGG_2021_Human"],
            gsea_library="KEGG_2021_Human",
        )
        result = Runner(dataset, config).run_all()
        from app.core.report import build_report

        report_path = build_report(result, dataset, config)
        report["run_dir"] = str(result.project_dir)
        report["figures"] = len(result.figures)
        report["tables"] = len(result.tables)
        report["report"] = str(report_path)
        report["report_bytes"] = Path(report_path).stat().st_size
        report["errors"] = result.errors
    except Exception as exc:  # 自检失败也要留证据
        import traceback

        report["fatal"] = f"{type(exc).__name__}: {exc}"
        report["traceback"] = traceback.format_exc()

    Path(report_file).write_text(json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
    return 0 if "fatal" not in report else 1


def main() -> int:
    if len(sys.argv) >= 3 and sys.argv[1] == "--selftest":
        return _selftest(sys.argv[2])

    # Qt 6 已默认启用高 DPI 缩放与高清位图；AA_EnableHighDpiScaling / AA_UseHighDpiPixmaps
    # 在 Qt 6 中被标记废弃且不再生效，故不再设置。
    app = QApplication(sys.argv)
    app.setApplicationName("OmicsOne")
    app.setApplicationDisplayName("OmicsOne 一键组学分析平台")
    app.setStyle("Fusion")
    app.setFont(QFont("Microsoft YaHei UI", 10))
    app.setStyleSheet(APP_QSS)

    icon_file = logo_path()
    if icon_file:
        app.setWindowIcon(QIcon(str(icon_file)))

    sys.excepthook = _excepthook

    window = MainWindow()
    window.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
