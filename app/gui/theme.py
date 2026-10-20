"""暗色科技风主题：配色、字体与全局 QSS。"""

from __future__ import annotations

from pathlib import Path

# 主色板（与报告、图表保持一致）
BG_DEEP = "#070E16"
BG_BASE = "#0B1622"
BG_PANEL = "#101E2B"
BG_PANEL_ALT = "#132434"
BORDER = "#1E3346"
BORDER_SOFT = "#1A2B3B"
TEXT = "#E6F1F6"
TEXT_DIM = "#8FA6B8"
TEXT_MUTE = "#5F7787"
ACCENT = "#2FD8C6"
ACCENT_BLUE = "#3FA9F5"
ACCENT_DEEP = "#1C8F9C"
DANGER = "#F2707A"
WARN = "#F4C77B"
OK = "#5FD68C"

FONT_FAMILY = "Microsoft YaHei UI, Microsoft YaHei, Segoe UI, sans-serif"

APP_QSS = f"""
* {{
    font-family: {FONT_FAMILY};
    outline: none;
}}
QWidget#Root {{
    background: {BG_BASE};
}}
QWidget#TopBar {{
    background: {BG_PANEL};
    border-bottom: 1px solid {BORDER_SOFT};
}}
QWidget#SideBar {{
    background: {BG_PANEL};
    border-right: 1px solid {BORDER_SOFT};
}}
QWidget#Content {{
    background: {BG_BASE};
}}
QLabel#Brand {{
    color: {TEXT};
    font-size: 19px;
    font-weight: 600;
    letter-spacing: 0.3px;
}}
QLabel#BrandMark {{
    border-radius: 9px;
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
        stop:0 {ACCENT_BLUE}, stop:1 {ACCENT});
    color: #05161C;
    font-size: 16px;
    font-weight: 700;
}}
QFrame#Card {{
    background: {BG_PANEL};
    border: 1px solid {BORDER};
    border-radius: 12px;
}}
QFrame#CardHeader {{
    background: transparent;
    border: none;
    border-bottom: 1px solid {BORDER_SOFT};
}}
QLabel#CardTitle {{
    color: {TEXT};
    font-size: 15px;
    font-weight: 600;
}}
QLabel#SectionLabel {{
    color: {TEXT_DIM};
    font-size: 12px;
}}
QLabel#Hint {{
    color: {TEXT_MUTE};
    font-size: 12px;
}}
QLabel#Value {{
    color: {TEXT};
    font-size: 13px;
}}
QLabel#BigValue {{
    color: {TEXT};
    font-size: 26px;
    font-weight: 600;
}}
QPushButton {{
    background: {BG_PANEL_ALT};
    color: {TEXT};
    border: 1px solid {BORDER};
    border-radius: 7px;
    padding: 7px 16px;
    font-size: 13px;
}}
QPushButton:hover {{
    border-color: {ACCENT_DEEP};
    background: #172E3F;
}}
QPushButton:pressed {{
    background: #0E2231;
}}
QPushButton:disabled {{
    color: {TEXT_MUTE};
    border-color: {BORDER_SOFT};
    background: #0D1A26;
}}
QPushButton#Primary {{
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
        stop:0 {ACCENT_BLUE}, stop:0.55 #33C7DE, stop:1 {ACCENT});
    color: #04161B;
    border: none;
    border-radius: 9px;
    font-size: 19px;
    font-weight: 700;
    letter-spacing: 1px;
    padding: 14px 24px;
}}
QPushButton#Primary:hover {{
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
        stop:0 #55BBFF, stop:0.55 #45DCEC, stop:1 #52EDD9);
}}
QPushButton#Primary:disabled {{
    background: #17303F;
    color: #4E6B7C;
}}
QPushButton#Danger {{
    background: transparent;
    color: {DANGER};
    border: 1px solid #54303A;
}}
QPushButton#Danger:hover {{
    background: #2A1A22;
    border-color: {DANGER};
}}
QPushButton#Ghost {{
    background: transparent;
    border: 1px solid {BORDER};
    color: {TEXT_DIM};
}}
QPushButton#Ghost:hover {{
    color: {TEXT};
    border-color: {ACCENT_DEEP};
}}
QPushButton#Nav {{
    background: transparent;
    border: none;
    border-left: 3px solid transparent;
    border-radius: 0;
    color: {TEXT_DIM};
    text-align: left;
    padding: 13px 18px;
    font-size: 14px;
}}
QPushButton#Nav:hover {{
    color: {TEXT};
    background: #0F1E2B;
}}
QPushButton#Nav:checked {{
    color: {TEXT};
    background: #12293A;
    border-left: 3px solid {ACCENT};
    font-weight: 600;
}}
QPushButton#SegOn {{
    background: #DCE6EC;
    color: #0B1622;
    border: none;
    border-radius: 15px;
    padding: 6px 20px;
    font-size: 13px;
    font-weight: 600;
}}
QPushButton#SegOff {{
    background: transparent;
    color: {TEXT_DIM};
    border: none;
    border-radius: 15px;
    padding: 6px 18px;
    font-size: 13px;
}}
QPushButton#SegOff:hover {{
    color: {TEXT};
}}
QLineEdit, QComboBox, QSpinBox, QDoubleSpinBox {{
    background: #0C1926;
    border: 1px solid {BORDER};
    border-radius: 7px;
    color: {TEXT};
    padding: 7px 11px;
    font-size: 13px;
    selection-background-color: {ACCENT_DEEP};
}}
QLineEdit:focus, QComboBox:focus, QSpinBox:focus, QDoubleSpinBox:focus {{
    border-color: {ACCENT};
}}
QLineEdit:read-only {{
    color: {TEXT_DIM};
}}
QComboBox::drop-down {{
    border: none;
    width: 22px;
}}
QComboBox QAbstractItemView {{
    background: {BG_PANEL_ALT};
    border: 1px solid {BORDER};
    color: {TEXT};
    selection-background-color: {ACCENT_DEEP};
    padding: 4px;
    outline: none;
}}
QCheckBox, QRadioButton {{
    color: {TEXT_DIM};
    spacing: 7px;
    font-size: 13px;
}}
QCheckBox:hover, QRadioButton:hover {{
    color: {TEXT};
}}
QCheckBox::indicator, QRadioButton::indicator {{
    width: 15px;
    height: 15px;
    border: 1px solid #35505F;
    background: #0C1926;
}}
QCheckBox::indicator {{
    border-radius: 4px;
}}
QRadioButton::indicator {{
    border-radius: 8px;
}}
QCheckBox::indicator:checked {{
    background: {ACCENT};
    border-color: {ACCENT};
    image: none;
}}
QRadioButton::indicator:checked {{
    background: {ACCENT};
    border-color: {ACCENT};
}}
QCheckBox::indicator:disabled, QRadioButton::indicator:disabled {{
    border-color: #24333F;
    background: #0A141D;
}}
QProgressBar {{
    background: #0B1927;
    border: none;
    border-radius: 4px;
    height: 9px;
    text-align: center;
    color: transparent;
}}
QProgressBar::chunk {{
    border-radius: 4px;
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
        stop:0 {ACCENT_BLUE}, stop:1 {ACCENT});
}}
QTextEdit#LogView {{
    background: #08131C;
    border: 1px solid {BORDER_SOFT};
    border-radius: 8px;
    color: #C6DAE4;
    font-family: Consolas, "Cascadia Mono", monospace;
    font-size: 12px;
    padding: 8px;
}}
QScrollBar:vertical {{
    background: transparent;
    width: 9px;
    margin: 2px;
}}
QScrollBar::handle:vertical {{
    background: #23414F;
    border-radius: 4px;
    min-height: 26px;
}}
QScrollBar::handle:vertical:hover {{
    background: #2F5A6B;
}}
QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {{
    height: 0;
}}
QScrollBar:horizontal {{
    background: transparent;
    height: 9px;
    margin: 2px;
}}
QScrollBar::handle:horizontal {{
    background: #23414F;
    border-radius: 4px;
    min-width: 26px;
}}
QScrollBar::add-line:horizontal, QScrollBar::sub-line:horizontal {{
    width: 0;
}}
QScrollArea {{
    background: transparent;
    border: none;
}}
QScrollArea > QWidget > QWidget {{
    background: transparent;
}}
QToolTip {{
    background: {BG_PANEL_ALT};
    color: {TEXT};
    border: 1px solid {ACCENT_DEEP};
    border-radius: 6px;
    padding: 6px 9px;
}}
QSplitter::handle {{
    background: transparent;
}}
QTableWidget {{
    background: #0B1926;
    gridline-color: {BORDER_SOFT};
    color: {TEXT_DIM};
    border: 1px solid {BORDER_SOFT};
    border-radius: 8px;
    selection-background-color: #14303F;
}}
QHeaderView::section {{
    background: {BG_PANEL_ALT};
    color: {TEXT_DIM};
    border: none;
    border-bottom: 1px solid {BORDER};
    padding: 6px 8px;
    font-size: 12px;
}}
"""


def logo_path() -> Path | None:
    """打包后从 _MEIPASS 读，开发时从 app/assets 读。"""
    import sys

    base = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent.parent))
    for candidate in (base / "assets" / "logo.png", Path(__file__).resolve().parent.parent / "assets" / "logo.png"):
        if candidate.exists():
            return candidate
    return None
