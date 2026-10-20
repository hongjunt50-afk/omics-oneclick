"""界面截图工具：逐页截图，供视觉验收与文档使用。"""

from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from PySide6.QtCore import QTimer
from PySide6.QtGui import QFont
from PySide6.QtWidgets import QApplication

from app.core.data_loader import load_dataset
from app.gui import main_window as mw
from app.gui.main_window import MainWindow
from app.gui.theme import APP_QSS


def capture(out_dir: Path, data: str, groups: str | None) -> None:
    out_dir.mkdir(parents=True, exist_ok=True)
    app = QApplication(sys.argv)
    app.setStyle("Fusion")
    app.setFont(QFont("Microsoft YaHei UI", 10))
    app.setStyleSheet(APP_QSS)

    window = MainWindow()
    window.resize(1440, 900)
    window.show()

    dataset = load_dataset(data, group_file=groups)
    window.output_edit.setText(str(Path("OmicsOne_Results").resolve()))
    window.output_label.setText(str(Path("OmicsOne_Results").resolve()))
    window._on_dataset_ready(dataset, data, True)

    pages = {
        "01_一键全流程": mw.PAGE_ONE_CLICK,
        "02_数据导入": mw.PAGE_IMPORT,
        "03_输出设置": mw.PAGE_OUTPUT,
        "04_运行日志": mw.PAGE_LOG,
    }

    def shoot() -> None:
        for name, index in pages.items():
            window._go(index)
            app.processEvents()
            window.grab().save(str(out_dir / f"{name}.png"))
        app.quit()

    QTimer.singleShot(900, shoot)
    app.exec()
    print(f"截图输出到 {out_dir}")


if __name__ == "__main__":
    base = Path(__file__).resolve().parent.parent
    capture(
        base / "shots",
        str(base / "examples" / "demo_expression_matrix.csv"),
        str(base / "examples" / "demo_phenotype.csv"),
    )
