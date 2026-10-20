"""可复用界面组件：卡片、参数行、统计徽章、日志面板、进度条。"""

from __future__ import annotations

from PySide6.QtCore import Qt, Signal
from PySide6.QtWidgets import (
    QFrame,
    QHBoxLayout,
    QLabel,
    QPushButton,
    QSizePolicy,
    QTextEdit,
    QVBoxLayout,
    QWidget,
)

from . import theme


class Card(QFrame):
    """带标题栏的圆角卡片，内容放进 body 布局。"""

    def __init__(self, title: str = "", *, subtitle: str = "", parent: QWidget | None = None):
        super().__init__(parent)
        self.setObjectName("Card")
        outer = QVBoxLayout(self)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        self.header = QFrame()
        self.header.setObjectName("CardHeader")
        header_layout = QHBoxLayout(self.header)
        header_layout.setContentsMargins(20, 15, 20, 13)
        header_layout.setSpacing(10)

        self.title_label = QLabel(title)
        self.title_label.setObjectName("CardTitle")
        header_layout.addWidget(self.title_label)
        if subtitle:
            sub = QLabel(subtitle)
            sub.setObjectName("Hint")
            header_layout.addWidget(sub)
        header_layout.addStretch(1)
        self.header_layout = header_layout
        outer.addWidget(self.header)

        self.body = QWidget()
        self.body_layout = QVBoxLayout(self.body)
        self.body_layout.setContentsMargins(20, 18, 20, 18)
        self.body_layout.setSpacing(14)
        outer.addWidget(self.body, 1)

    def add_header_widget(self, widget: QWidget) -> None:
        self.header_layout.addWidget(widget)

    def add(self, widget: QWidget, stretch: int = 0) -> QWidget:
        self.body_layout.addWidget(widget, stretch)
        return widget

    def add_row(self, *widgets: QWidget, spacing: int = 10) -> QHBoxLayout:
        row = QHBoxLayout()
        row.setSpacing(spacing)
        for widget in widgets:
            if isinstance(widget, int):
                row.addStretch(widget)
            else:
                row.addWidget(widget)
        self.body_layout.addLayout(row)
        return row


class StatChip(QWidget):
    """顶部或卡片里的小指标块。"""

    def __init__(self, label: str, value: str = "—", parent: QWidget | None = None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(2)
        self.value_label = QLabel(value)
        self.value_label.setObjectName("BigValue")
        self.caption = QLabel(label)
        self.caption.setObjectName("SectionLabel")
        layout.addWidget(self.value_label)
        layout.addWidget(self.caption)

    def set_value(self, value: str) -> None:
        self.value_label.setText(value)


class Badge(QLabel):
    """状态徽章：ok / warn / err / idle。"""

    COLORS = {
        "ok": (theme.ACCENT, "#08231F"),
        "warn": (theme.WARN, "#2A2312"),
        "err": (theme.DANGER, "#2A161A"),
        "idle": (theme.TEXT_MUTE, "#12202C"),
    }

    def __init__(self, text: str = "", kind: str = "idle", parent: QWidget | None = None):
        super().__init__(text, parent)
        self.setAlignment(Qt.AlignCenter)
        self.set_state(kind, text)

    def set_state(self, kind: str, text: str | None = None) -> None:
        fg, bg = self.COLORS.get(kind, self.COLORS["idle"])
        if text is not None:
            self.setText(text)
        self.setStyleSheet(
            f"color:{fg}; background:{bg}; border:1px solid {fg}33;"
            "border-radius:6px; padding:3px 10px; font-size:12px;"
        )


class LogPanel(QWidget):
    """带级别着色的日志区，支持清空与导出。"""

    cleared = Signal()

    def __init__(self, parent: QWidget | None = None):
        super().__init__(parent)
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(10)

        self.view = QTextEdit()
        self.view.setObjectName("LogView")
        self.view.setReadOnly(True)
        self.view.setMinimumHeight(150)
        layout.addWidget(self.view, 1)

        buttons = QHBoxLayout()
        buttons.setSpacing(8)
        buttons.addStretch(1)
        self.clear_btn = QPushButton("清空日志")
        self.clear_btn.setObjectName("Ghost")
        self.export_btn = QPushButton("导出日志")
        self.export_btn.setObjectName("Ghost")
        buttons.addWidget(self.clear_btn)
        buttons.addWidget(self.export_btn)
        layout.addLayout(buttons)

        self.clear_btn.clicked.connect(self._clear)

    def _clear(self) -> None:
        self.view.clear()
        self.cleared.emit()

    def append(self, level: str, message: str, timestamp: str = "") -> None:
        color = {
            "ERROR": theme.DANGER,
            "WARNING": theme.WARN,
            "SUCCESS": theme.OK,
            "INFO": theme.ACCENT,
        }.get(level.upper(), theme.TEXT_DIM)
        stamp = f'<span style="color:{theme.TEXT_MUTE}">{timestamp}</span> ' if timestamp else ""
        safe = (
            message.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace("\n", "<br/>")
        )
        self.view.append(
            f'{stamp}<span style="color:{color}; font-weight:600">[{level.upper()}]</span> '
            f'<span style="color:#C6DAE4">{safe}</span>'
        )
        bar = self.view.verticalScrollBar()
        bar.setValue(bar.maximum())

    def text(self) -> str:
        return self.view.toPlainText()


class PathSelector(QWidget):
    """只读路径框 + 浏览按钮。"""

    def __init__(self, placeholder: str, button_text: str = "浏览…", parent: QWidget | None = None):
        super().__init__(parent)
        layout = QHBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(8)
        self.edit = QLabel(placeholder)
        self.edit.setObjectName("Value")
        self.edit.setStyleSheet(
            f"background:#0C1926; border:1px solid {theme.BORDER}; border-radius:7px;"
            f" padding:7px 11px; color:{theme.TEXT_DIM};"
        )
        self.edit.setMinimumHeight(34)
        self.edit.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.button = QPushButton(button_text)
        layout.addWidget(self.edit, 1)
        layout.addWidget(self.button)

    def set_text(self, text: str) -> None:
        self.edit.setText(text)
        self.edit.setToolTip(text)
        self.edit.setStyleSheet(
            f"background:#0C1926; border:1px solid {theme.BORDER}; border-radius:7px;"
            f" padding:7px 11px; color:{theme.TEXT};"
        )
