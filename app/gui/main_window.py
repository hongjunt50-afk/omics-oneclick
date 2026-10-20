"""主窗口：侧边导航 + 五个功能页，暗色科技风。"""

from __future__ import annotations

import time
from pathlib import Path

from PySide6.QtCore import Qt, QTimer, Signal
from PySide6.QtGui import QDesktopServices, QFont
from PySide6.QtCore import QUrl
from PySide6.QtWidgets import (
    QApplication,
    QButtonGroup,
    QCheckBox,
    QComboBox,
    QDoubleSpinBox,
    QFileDialog,
    QFrame,
    QGridLayout,
    QHBoxLayout,
    QHeaderView,
    QInputDialog,
    QLabel,
    QLineEdit,
    QMainWindow,
    QMessageBox,
    QProgressBar,
    QPushButton,
    QRadioButton,
    QScrollArea,
    QSizePolicy,
    QSpinBox,
    QStackedWidget,
    QTableWidget,
    QTableWidgetItem,
    QVBoxLayout,
    QWidget,
)

from .. import __version__
from ..core.data_loader import group_table_frame, save_group_table
from ..core.genesets import BUNDLED
from ..core.pipeline import AnalysisConfig
from . import theme
from .widgets import Badge, Card, LogPanel, StatChip
from .workers import AnalysisWorker, DownloadWorker, LoadWorker, open_folder

PAGE_ONE_CLICK = 0
PAGE_IMPORT = 1
PAGE_OUTPUT = 2
PAGE_LOG = 3

# (键, 显示名, 所属流水线步骤, 说明)
ANALYSIS_ITEMS = [
    ("qc", "质控 QC", "qc", "样本分布、相关性、离群样本检测"),
    ("norm", "标准化", "qc", "分位数 / 中位数标准化"),
    ("deg", "差异分析 DEG", "deg", "limma 风格线性模型 + 经验贝叶斯"),
    ("pca", "PCA 降维", "pca", "主成分分析与样本分离度"),
    ("go", "GO 富集", "enrichment", "超几何检验 + BH 校正"),
    ("kegg", "KEGG 富集", "enrichment", "通路过表达分析"),
    ("volcano", "火山图", "deg", "差异倍数与显著性分布"),
    ("heatmap", "聚类热图", "deg", "Top 差异基因表达模式"),
    ("wgcna", "WGCNA", "wgcna", "共表达模块与性状关联"),
    ("gsea", "GSEA", "gsea", "预排序基因集富集分析"),
]

STEP_LABELS = {
    "qc": "质控与标准化",
    "deg": "差异基因分析",
    "pca": "主成分分析",
    "enrichment": "富集分析",
    "gsea": "GSEA 分析",
    "wgcna": "WGCNA 分析",
}


class MainWindow(QMainWindow):
    def __init__(self) -> None:
        super().__init__()
        self.setWindowTitle(f"OmicsOne 一键组学分析平台  v{__version__}")
        self.resize(1320, 880)
        self.setMinimumSize(1120, 760)

        self.dataset = None
        self.dataset_path: Path | None = None
        self.group_file: str | None = None
        self.worker: AnalysisWorker | None = None
        self.load_worker: LoadWorker | None = None
        self.download_worker: DownloadWorker | None = None
        self.last_result = None
        self.last_report: str | None = None
        self._started_at = 0.0
        self._last_fraction = 0.0
        self._current_task = "等待开始"
        self._paused_at = 0.0

        self._eta_timer = QTimer(self)
        self._eta_timer.setInterval(1000)
        self._eta_timer.timeout.connect(self._tick_eta)

        self._build_ui()
        self._apply_state()
        self._log("INFO", "OmicsOne 已启动，请先导入数据集")

    # ==================== 界面搭建 ====================
    def _build_ui(self) -> None:
        root = QWidget()
        root.setObjectName("Root")
        self.setCentralWidget(root)
        outer = QVBoxLayout(root)
        outer.setContentsMargins(0, 0, 0, 0)
        outer.setSpacing(0)

        outer.addWidget(self._build_topbar())

        body = QHBoxLayout()
        body.setContentsMargins(0, 0, 0, 0)
        body.setSpacing(0)
        body.addWidget(self._build_sidebar())

        self.stack = QStackedWidget()
        self.stack.setObjectName("Content")
        self.stack.addWidget(self._wrap_scroll(self._build_page_one_click()))
        self.stack.addWidget(self._wrap_scroll(self._build_page_import()))
        self.stack.addWidget(self._wrap_scroll(self._build_page_output()))
        self.stack.addWidget(self._build_page_log())
        body.addWidget(self.stack, 1)

        holder = QWidget()
        holder.setObjectName("Content")
        holder.setLayout(body)
        outer.addWidget(holder, 1)

        self._sync_checks_from_types()

    def _wrap_scroll(self, widget: QWidget) -> QScrollArea:
        area = QScrollArea()
        area.setWidgetResizable(True)
        area.setWidget(widget)
        area.setFrameShape(QFrame.NoFrame)
        area.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        return area

    def _build_topbar(self) -> QWidget:
        bar = QWidget()
        bar.setObjectName("TopBar")
        bar.setFixedHeight(66)
        layout = QHBoxLayout(bar)
        layout.setContentsMargins(20, 0, 22, 0)
        layout.setSpacing(12)

        mark = QLabel("◈")
        mark.setObjectName("BrandMark")
        mark.setFixedSize(30, 30)
        mark.setAlignment(Qt.AlignCenter)
        layout.addWidget(mark)

        title = QLabel("OmicsOne  一键组学分析平台")
        title.setObjectName("Brand")
        layout.addWidget(title)

        layout.addStretch(1)

        self.project_combo = QComboBox()
        self.project_combo.setMinimumWidth(300)
        self.project_combo.addItem("当前项目：未命名项目")
        layout.addWidget(self.project_combo)

        layout.addSpacing(8)

        self.top_badge = Badge("等待数据", "idle")
        layout.addWidget(self.top_badge)

        seg_file = QPushButton("文件")
        seg_file.setObjectName("SegOn")
        seg_file.setCursor(Qt.PointingHandCursor)
        seg_file.clicked.connect(self._pick_file)
        layout.addWidget(seg_file)

        settings_btn = QPushButton("设置")
        settings_btn.setObjectName("Ghost")
        settings_btn.clicked.connect(lambda: self._go(PAGE_OUTPUT))
        layout.addWidget(settings_btn)

        help_btn = QPushButton("帮助")
        help_btn.setObjectName("Ghost")
        help_btn.clicked.connect(self._show_help)
        layout.addWidget(help_btn)

        bell = QLabel("◔")
        bell.setToolTip(f"OmicsOne v{__version__}")
        bell.setStyleSheet(f"color:{theme.TEXT_DIM}; font-size:17px; padding:0 4px;")
        layout.addWidget(bell)

        avatar = QLabel("O1")
        avatar.setFixedSize(30, 30)
        avatar.setAlignment(Qt.AlignCenter)
        avatar.setStyleSheet(
            f"background:{theme.ACCENT_DEEP}; color:#04161B; border-radius:15px;"
            "font-size:11px; font-weight:700;"
        )
        layout.addWidget(avatar)
        return bar

    def _build_sidebar(self) -> QWidget:
        side = QWidget()
        side.setObjectName("SideBar")
        side.setFixedWidth(196)
        layout = QVBoxLayout(side)
        layout.setContentsMargins(0, 14, 0, 14)
        layout.setSpacing(2)

        self.nav_group = QButtonGroup(self)
        self.nav_group.setExclusive(True)
        entries = [
            ("◈  数据导入", PAGE_IMPORT),
            ("▤  一键全流程分析", PAGE_ONE_CLICK),
            ("⌸  结果输出设置", PAGE_OUTPUT),
            ("≡  运行日志", PAGE_LOG),
        ]
        for text, index in entries:
            btn = QPushButton(text)
            btn.setObjectName("Nav")
            btn.setCheckable(True)
            btn.setCursor(Qt.PointingHandCursor)
            btn.clicked.connect(lambda _=False, i=index: self._go(i))
            self.nav_group.addButton(btn, index)
            layout.addWidget(btn)

        layout.addStretch(1)
        tip = QLabel("全部计算在本机完成\n数据不上传、不联网")
        tip.setObjectName("Hint")
        tip.setAlignment(Qt.AlignCenter)
        tip.setWordWrap(True)
        layout.addWidget(tip)
        self.nav_group.button(PAGE_ONE_CLICK).setChecked(True)
        return side

    def _page(self) -> QWidget:
        page = QWidget()
        page.setObjectName("Content")
        layout = QVBoxLayout(page)
        layout.setContentsMargins(26, 22, 26, 26)
        layout.setSpacing(18)
        return page

    # ---------- 页 1：一键全流程 ----------
    def _build_page_one_click(self) -> QWidget:
        page = self._page()
        layout = page.layout()

        header = QHBoxLayout()
        title = QLabel("一键组学全流程分析")
        title.setStyleSheet(f"color:{theme.TEXT}; font-size:21px; font-weight:600;")
        header.addWidget(title)
        header.addStretch(1)
        out_caption = QLabel("结果输出文件夹：")
        out_caption.setObjectName("SectionLabel")
        header.addWidget(out_caption)
        self.output_label = QLabel("未设置")
        self.output_label.setObjectName("Value")
        header.addWidget(self.output_label)
        pick = QPushButton("更改")
        pick.setObjectName("Ghost")
        pick.clicked.connect(self._choose_output_dir)
        header.addWidget(pick)
        layout.addLayout(header)

        # --- 数据集导入卡片 ---
        import_card = Card("数据集导入")
        self.import_badge = Badge("未导入", "idle")
        import_card.add_header_widget(self.import_badge)

        btns = QHBoxLayout()
        btns.setSpacing(10)
        self.upload_btn = QPushButton("本地上传")
        self.upload_btn.setObjectName("Primary")
        self.upload_btn.setStyleSheet(
            f"background:{theme.ACCENT}; color:#04231F; border:none; border-radius:7px;"
            "padding:9px 22px; font-size:13px; font-weight:600;"
        )
        self.url_btn = QPushButton("链接下载")
        self.folder_btn = QPushButton("选择文件夹")
        self.group_btn = QPushButton("分组表（可选）")
        self.group_btn.setObjectName("Ghost")
        self.auto_group_btn = QPushButton("自动生成分组表")
        self.auto_group_btn.setObjectName("Ghost")
        self.auto_group_btn.setToolTip("按样本名自动推断分组，导出 sample,group 两列 CSV 并直接载入")
        for button in (self.upload_btn, self.url_btn, self.folder_btn, self.auto_group_btn):
            button.setCursor(Qt.PointingHandCursor)
        btns.addWidget(self.upload_btn)
        btns.addWidget(self.url_btn)
        btns.addWidget(self.folder_btn)
        btns.addWidget(self.group_btn)
        btns.addWidget(self.auto_group_btn)
        btns.addStretch(1)
        import_card.body_layout.addLayout(btns)

        file_row = QHBoxLayout()
        file_row.setSpacing(10)
        self.file_label = QLabel("尚未选择数据集文件")
        self.file_label.setObjectName("Value")
        self.file_label.setStyleSheet(
            f"background:#0C1926; border:1px solid {theme.BORDER}; border-radius:7px;"
            f" padding:8px 12px; color:{theme.TEXT_MUTE}; min-height:20px;"
        )
        self.file_label.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        file_row.addWidget(self.file_label, 1)
        self.file_badge = Badge("待加载", "idle")
        file_row.addWidget(self.file_badge)
        self.preview_btn = QPushButton("预览数据")
        self.preview_btn.setObjectName("Ghost")
        self.preview_btn.setEnabled(False)
        file_row.addWidget(self.preview_btn)
        import_card.body_layout.addLayout(file_row)

        self.dataset_summary = QLabel("导入后这里显示基因数、样本数与分组信息")
        self.dataset_summary.setObjectName("Hint")
        self.dataset_summary.setWordWrap(True)
        import_card.body_layout.addWidget(self.dataset_summary)
        layout.addWidget(import_card)

        # --- 分析配置卡片 ---
        config_card = Card("分析配置")

        type_row = QHBoxLayout()
        type_row.setSpacing(14)
        type_caption = QLabel("组学类型")
        type_caption.setObjectName("SectionLabel")
        type_row.addWidget(type_caption)
        self.type_group = QButtonGroup(self)
        for index, name in enumerate(("转录组", "蛋白组", "代谢组")):
            radio = QRadioButton(name)
            if index == 0:
                radio.setChecked(True)
            radio.toggled.connect(self._on_type_changed)
            self.type_group.addButton(radio, index)
            type_row.addWidget(radio)
        type_row.addSpacing(18)

        cmp_caption = QLabel("比较组")
        cmp_caption.setObjectName("SectionLabel")
        type_row.addWidget(cmp_caption)
        self.group_a_combo = QComboBox()
        self.group_b_combo = QComboBox()
        self.group_a_combo.setMinimumWidth(120)
        self.group_b_combo.setMinimumWidth(120)
        type_row.addWidget(self.group_a_combo)
        type_row.addWidget(QLabel("vs"))
        type_row.addWidget(self.group_b_combo)
        type_row.addStretch(1)
        config_card.body_layout.addLayout(type_row)

        self.checkboxes: dict[str, QCheckBox] = {}
        grid = QGridLayout()
        grid.setHorizontalSpacing(26)
        grid.setVerticalSpacing(10)
        for i, (key, label, _step, tip) in enumerate(ANALYSIS_ITEMS):
            box = QCheckBox(label)
            box.setToolTip(tip)
            box.setChecked(True)
            box.toggled.connect(self._sync_run_button)
            self.checkboxes[key] = box
            grid.addWidget(box, i // 5, i % 5)
        config_card.body_layout.addLayout(grid)

        option_row = QHBoxLayout()
        option_row.setSpacing(16)
        padj_caption = QLabel("padj <")
        padj_caption.setObjectName("SectionLabel")
        self.padj_spin = QDoubleSpinBox()
        self.padj_spin.setDecimals(3)
        self.padj_spin.setSingleStep(0.01)
        self.padj_spin.setRange(0.001, 0.5)
        self.padj_spin.setValue(0.05)
        self.padj_spin.setFixedWidth(90)

        lfc_caption = QLabel("|log2FC| ≥")
        lfc_caption.setObjectName("SectionLabel")
        self.lfc_spin = QDoubleSpinBox()
        self.lfc_spin.setDecimals(2)
        self.lfc_spin.setSingleStep(0.1)
        self.lfc_spin.setRange(0.0, 5.0)
        self.lfc_spin.setValue(1.0)
        self.lfc_spin.setFixedWidth(90)

        norm_caption = QLabel("标准化")
        norm_caption.setObjectName("SectionLabel")
        self.norm_combo = QComboBox()
        self.norm_combo.addItem("分位数标准化", "quantile")
        self.norm_combo.addItem("中位数中心化", "median")
        self.norm_combo.addItem("不做标准化", "none")
        self.norm_combo.setFixedWidth(140)

        perm_caption = QLabel("GSEA 置换")
        perm_caption.setObjectName("SectionLabel")
        self.perm_spin = QSpinBox()
        self.perm_spin.setRange(100, 10000)
        self.perm_spin.setSingleStep(100)
        self.perm_spin.setValue(1000)
        self.perm_spin.setFixedWidth(90)

        for widget in (
            padj_caption, self.padj_spin, lfc_caption, self.lfc_spin,
            norm_caption, self.norm_combo, perm_caption, self.perm_spin,
        ):
            option_row.addWidget(widget)
        option_row.addStretch(1)
        config_card.body_layout.addLayout(option_row)

        action_row = QHBoxLayout()
        action_row.setSpacing(12)
        self.run_btn = QPushButton("🚀  一键开始分析")
        self.run_btn.setObjectName("Primary")
        self.run_btn.setMinimumHeight(58)
        self.run_btn.setCursor(Qt.PointingHandCursor)
        self.run_btn.clicked.connect(self._run_one_click)
        action_row.addWidget(self.run_btn, 3)

        self.pause_btn = QPushButton("暂停")
        self.pause_btn.setEnabled(False)
        self.pause_btn.clicked.connect(self._toggle_pause)
        self.stop_btn = QPushButton("终止")
        self.stop_btn.setObjectName("Danger")
        self.stop_btn.setEnabled(False)
        self.stop_btn.clicked.connect(self._stop_run)
        action_row.addWidget(self.pause_btn, 1)
        action_row.addWidget(self.stop_btn, 1)
        config_card.body_layout.addLayout(action_row)
        layout.addWidget(config_card)

        # --- 进度卡片 ---
        progress_card = Card("运行进度")
        row = QHBoxLayout()
        self.task_label = QLabel("当前任务：等待开始")
        self.task_label.setObjectName("Value")
        row.addWidget(self.task_label)
        row.addStretch(1)
        self.percent_label = QLabel("0%")
        self.percent_label.setStyleSheet(f"color:{theme.ACCENT}; font-size:16px; font-weight:600;")
        row.addWidget(self.percent_label)
        row.addSpacing(18)
        self.eta_label = QLabel("ETA --:--")
        self.eta_label.setObjectName("Hint")
        row.addWidget(self.eta_label)
        progress_card.body_layout.addLayout(row)

        self.progress_bar = QProgressBar()
        self.progress_bar.setRange(0, 1000)
        self.progress_bar.setValue(0)
        self.progress_bar.setTextVisible(False)
        self.progress_bar.setFixedHeight(9)
        progress_card.body_layout.addWidget(self.progress_bar)

        result_row = QHBoxLayout()
        self.result_label = QLabel("分析完成后可一键打开报告或结果目录")
        self.result_label.setObjectName("Hint")
        result_row.addWidget(self.result_label)
        result_row.addStretch(1)
        self.report_btn = QPushButton("打开报告")
        self.report_btn.setEnabled(False)
        self.report_btn.clicked.connect(self._open_report)
        self.folder_result_btn = QPushButton("打开结果目录")
        self.folder_result_btn.setEnabled(False)
        self.folder_result_btn.clicked.connect(self._open_result_folder)
        result_row.addWidget(self.report_btn)
        result_row.addWidget(self.folder_result_btn)
        progress_card.body_layout.addLayout(result_row)
        layout.addWidget(progress_card)

        # --- 日志卡片 ---
        log_card = Card("运行日志")
        self.one_click_log = LogPanel()
        self.one_click_log.setMinimumHeight(170)
        self.one_click_log.export_btn.clicked.connect(self._export_log)
        log_card.add(self.one_click_log)
        layout.addWidget(log_card)
        layout.addStretch(1)

        self.upload_btn.clicked.connect(self._pick_file)
        self.url_btn.clicked.connect(self._download_dialog)
        self.folder_btn.clicked.connect(self._pick_folder)
        self.group_btn.clicked.connect(self._pick_group_file)
        self.auto_group_btn.clicked.connect(self._auto_group_table)
        self.preview_btn.clicked.connect(lambda: self._go(PAGE_IMPORT))
        return page

    # ---------- 页 2：数据导入 ----------
    def _build_page_import(self) -> QWidget:
        page = self._page()
        layout = page.layout()

        title = QLabel("数据导入与检查")
        title.setStyleSheet(f"color:{theme.TEXT}; font-size:21px; font-weight:600;")
        layout.addWidget(title)

        source_card = Card("数据来源")
        row = QHBoxLayout()
        row.setSpacing(10)
        self.import_local_btn = QPushButton("本地上传文件")
        self.import_url_btn = QPushButton("通过链接下载")
        self.import_folder_btn = QPushButton("选择文件夹")
        self.import_auto_group_btn = QPushButton("自动生成分组表")
        self.import_auto_group_btn.setObjectName("Ghost")
        self.import_auto_group_btn.setToolTip("按样本名自动推断分组，导出 sample,group 两列 CSV 并直接载入")
        for button in (self.import_local_btn, self.import_url_btn, self.import_folder_btn,
                       self.import_auto_group_btn):
            row.addWidget(button)
        row.addStretch(1)
        source_card.body_layout.addLayout(row)
        self.import_path_label = QLabel("尚未导入数据")
        self.import_path_label.setObjectName("Hint")
        self.import_path_label.setWordWrap(True)
        source_card.body_layout.addWidget(self.import_path_label)
        layout.addWidget(source_card)

        # 手动填写 URL
        link_card = Card("链接下载")
        link_card.add_header_widget(Badge("支持直链 / 压缩包", "idle"))
        url_row = QHBoxLayout()
        self.url_edit = QLineEdit()
        self.url_edit.setPlaceholderText("粘贴数据直链，例如 https://example.com/GSE12345_matrix.csv.gz")
        self.url_go_btn = QPushButton("开始下载")
        self.url_go_btn.setObjectName("Ghost")
        url_row.addWidget(self.url_edit, 1)
        url_row.addWidget(self.url_go_btn)
        link_card.body_layout.addLayout(url_row)
        self.url_hint = QLabel("支持 csv / tsv / txt / tab / xlsx / xls / csv.gz，以及包含这些文件的 zip / tar.gz；下载后自动解压并识别。")
        self.url_hint.setObjectName("Hint")
        self.url_hint.setWordWrap(True)
        link_card.body_layout.addWidget(self.url_hint)
        self.url_progress = QProgressBar()
        self.url_progress.setRange(0, 1000)
        self.url_progress.setValue(0)
        self.url_progress.setTextVisible(False)
        self.url_progress.setFixedHeight(6)
        link_card.body_layout.addWidget(self.url_progress)
        layout.addWidget(link_card)

        preview_card = Card("数据预览")
        self.preview_info = QLabel("导入成功后这里显示前 10 行 × 8 列")
        self.preview_info.setObjectName("Hint")
        preview_card.add(self.preview_info)
        self.preview_table = QTableWidget()
        self.preview_table.setMinimumHeight(320)
        self.preview_table.horizontalHeader().setSectionResizeMode(QHeaderView.Interactive)
        self.preview_table.verticalHeader().setDefaultSectionSize(26)
        preview_card.add(self.preview_table, 1)
        layout.addWidget(preview_card, 1)

        self.import_local_btn.clicked.connect(self._pick_file)
        self.import_url_btn.clicked.connect(self._download_dialog)
        self.import_folder_btn.clicked.connect(self._pick_folder)
        self.import_auto_group_btn.clicked.connect(self._auto_group_table)
        self.url_go_btn.clicked.connect(self._download_from_edit)
        return page

    # ---------- 页 4：输出设置 ----------
    def _build_page_output(self) -> QWidget:
        page = self._page()
        layout = page.layout()

        title = QLabel("结果输出设置")
        title.setStyleSheet(f"color:{theme.TEXT}; font-size:21px; font-weight:600;")
        layout.addWidget(title)

        card = Card("输出位置")
        row = QHBoxLayout()
        row.setSpacing(10)
        self.output_edit = QLineEdit()
        self.output_edit.setPlaceholderText("选择结果保存文件夹")
        self.output_browse = QPushButton("浏览…")
        self.output_browse.clicked.connect(self._choose_output_dir)
        row.addWidget(self.output_edit, 1)
        row.addWidget(self.output_browse)
        card.body_layout.addLayout(row)

        name_row = QHBoxLayout()
        name_row.setSpacing(10)
        name_caption = QLabel("项目名")
        name_caption.setObjectName("SectionLabel")
        self.project_name_edit = QLineEdit("未命名项目")
        self.project_name_edit.setMaximumWidth(280)
        name_row.addWidget(name_caption)
        name_row.addWidget(self.project_name_edit)
        name_row.addStretch(1)
        card.body_layout.addLayout(name_row)
        layout.addWidget(card)

        struct_card = Card("输出目录结构")
        tree = QLabel(
            "OmicsOne_<项目名>_<时间戳>/\n"
            "├── 分析报告.html        自包含报告（图片内嵌，可直接发人）\n"
            "├── 分析概览.md          一页纸摘要\n"
            "├── summary.json        机器可读结果\n"
            "├── run_log.txt         完整运行日志\n"
            "├── figures/            全部图表 PNG（火山图、热图、PCA、富集、GSEA、WGCNA）\n"
            "└── tables/             全部结果表 CSV（UTF-8 BOM，Excel 双击可开）"
        )
        tree.setStyleSheet(
            f"color:{theme.TEXT_DIM}; font-family:Consolas, monospace; font-size:12px;"
            f"background:#0B1926; border:1px solid {theme.BORDER_SOFT}; border-radius:8px; padding:14px;"
        )
        struct_card.add(tree)
        layout.addWidget(struct_card)

        info_card = Card("数据集与分组")
        self.output_dataset_label = QLabel("尚未导入数据集")
        self.output_dataset_label.setObjectName("Hint")
        self.output_dataset_label.setWordWrap(True)
        info_card.add(self.output_dataset_label)
        self.output_group_btn = QPushButton("指定分组表…")
        self.output_group_btn.setObjectName("Ghost")
        self.output_group_btn.clicked.connect(self._pick_group_file)
        info_card.add(self.output_group_btn)
        layout.addWidget(info_card)
        layout.addStretch(1)
        return page

    # ---------- 页 5：日志 ----------
    def _build_page_log(self) -> QWidget:
        page = self._page()
        layout = page.layout()
        title = QLabel("运行日志")
        title.setStyleSheet(f"color:{theme.TEXT}; font-size:21px; font-weight:600;")
        layout.addWidget(title)

        card = Card("完整日志")
        self.log_panel = LogPanel()
        self.log_panel.setMinimumHeight(420)
        self.log_panel.export_btn.clicked.connect(self._export_log)
        card.add(self.log_panel)
        layout.addWidget(card, 1)
        return page

    # ==================== 导航与状态 ====================
    def _go(self, index: int) -> None:
        self.stack.setCurrentIndex(index)
        button = self.nav_group.button(index)
        if button is not None:
            button.setChecked(True)

    def _apply_state(self) -> None:
        has_data = self.dataset is not None
        busy = self.worker is not None and self.worker.isRunning()
        self.preview_btn.setEnabled(has_data)
        self.report_btn.setEnabled(bool(self.last_report))
        self.folder_result_btn.setEnabled(self.last_result is not None)
        self.stop_btn.setEnabled(busy)
        self.pause_btn.setEnabled(busy)
        for widget in (self.upload_btn, self.url_btn, self.folder_btn, self.import_local_btn,
                       self.import_url_btn, self.import_folder_btn, self.url_go_btn):
            widget.setEnabled(not busy)
        for widget in (self.auto_group_btn, self.import_auto_group_btn):
            widget.setEnabled(has_data and not busy)
        self._sync_run_button()

    def _sync_run_button(self) -> None:
        busy = self.worker is not None and self.worker.isRunning()
        selected = any(box.isChecked() for box in self.checkboxes.values())
        self.run_btn.setEnabled(self.dataset is not None and not busy and selected)
        if self.dataset is None:
            self.run_btn.setToolTip("请先导入数据集")
        elif not selected:
            self.run_btn.setToolTip("请至少勾选一项分析")
        else:
            self.run_btn.setToolTip("")

    def _sync_checks_from_types(self) -> None:
        pass

    def _on_type_changed(self) -> None:
        omics = self.type_group.checkedId()
        if omics == 2:  # 代谢组
            for key in ("go", "kegg", "gsea", "wgcna"):
                self.checkboxes[key].setChecked(False)
            self.checkboxes["volcano"].setChecked(True)
            self.lfc_spin.setValue(1.0)
            self._log("INFO", "已切换到代谢组：默认关闭基因符号依赖的分析（富集/GSEA/WGCNA）")
        else:
            for key in ("go", "kegg", "gsea", "wgcna", "volcano", "heatmap"):
                self.checkboxes[key].setChecked(True)
            self.lfc_spin.setValue(1.0 if omics == 0 else 1.0)

    # ==================== 数据导入 ====================
    def _pick_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self,
            "选择表达矩阵文件",
            "",
            "数据文件 (*.csv *.tsv *.txt *.tab *.xlsx *.xls *.gz);;所有文件 (*.*)",
        )
        if path:
            self._load_file(path)

    def _pick_folder(self) -> None:
        folder = QFileDialog.getExistingDirectory(self, "选择数据集文件夹")
        if not folder:
            return
        root = Path(folder)
        candidates = [
            p for p in sorted(root.iterdir())
            if p.is_file() and p.suffix.lower() in {".csv", ".tsv", ".txt", ".tab", ".xlsx", ".xls"}
        ]
        if not candidates:
            QMessageBox.warning(self, "未找到数据文件", "该文件夹下没有 csv/tsv/txt/xlsx 数据文件。")
            return
        if len(candidates) == 1:
            self._load_file(str(candidates[0]))
            return
        names = [p.name for p in candidates]
        choice, ok = QInputDialog.getItem(self, "选择数据文件", f"检测到 {len(names)} 个候选文件：", names, 0, False)
        if ok and choice:
            self._load_file(str(root / choice))

    def _pick_group_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "选择分组表（两列：样本名, 组名）", "", "表格 (*.csv *.tsv *.txt *.xlsx);;所有文件 (*.*)"
        )
        if not path:
            return
        self.group_file = path
        self._log("INFO", f"分组表：{Path(path).name}")
        if self.dataset_path:
            self._load_file(str(self.dataset_path), quiet=True)

    def _auto_group_table(self) -> None:
        """按样本名推断分组，落一张 sample,group 两列 CSV，并立即载入。"""
        if self.dataset is None:
            QMessageBox.information(self, "请先导入数据", "自动分组需要先导入表达矩阵。")
            return
        samples = [str(c) for c in self.dataset.expr.columns]
        frame, notes = group_table_frame(samples)
        if frame is None:
            QMessageBox.warning(
                self,
                "无法自动分组",
                "样本名里没有可识别的组别前缀（例如 control_1 / ipf_2）。\n\n"
                "请手动准备一张两列分组表（样本名, 组名），用「分组表（可选）」导入。",
            )
            return

        counts = frame["group"].value_counts()
        if len(counts) < 2 or int(counts.min()) < 2:
            QMessageBox.warning(
                self,
                "分组不可用",
                "自动推断出的分组不足两组，或某组样本少于 2 个：\n\n"
                + "、".join(f"{k}({v})" for k, v in counts.items())
                + "\n\n请手动准备一张两列分组表导入。",
            )
            return

        default = Path.home() / "OmicsOne_groups.csv"
        if self.dataset_path:
            default = Path(self.dataset_path).with_name(f"{Path(self.dataset_path).stem.split('.')[0]}_groups.csv")
        path, _ = QFileDialog.getSaveFileName(
            self, "保存分组表", str(default), "CSV (*.csv);;所有文件 (*.*)"
        )
        if not path:
            return

        try:
            saved = save_group_table(frame, path)
        except Exception as exc:
            QMessageBox.critical(self, "保存失败", str(exc))
            return

        self.group_file = str(saved)
        self._log("SUCCESS", f"已生成分组表：{saved.name}")
        self._log("INFO", "分组构成：" + "、".join(f"{k}({v})" for k, v in counts.items()))
        for note in notes:
            self._log("INFO", note)
        if self.dataset_path:
            self._load_file(str(self.dataset_path), quiet=True)
        else:
            QMessageBox.information(
                self, "分组表已保存", f"已写入 {saved}\n下次导入数据集时会自动带上这张分组表。"
            )

    def _load_file(self, path: str, *, quiet: bool = False) -> None:
        if self.worker is not None and self.worker.isRunning():
            QMessageBox.information(self, "请稍候", "分析正在运行，完成后再导入新数据。")
            return
        self.file_label.setText(Path(path).name)
        self.file_badge.set_state("warn", "加载中…")
        self.file_label.setStyleSheet(
            f"background:#0C1926; border:1px solid {theme.ACCENT}; border-radius:7px;"
            f" padding:8px 12px; color:{theme.TEXT}; min-height:20px;"
        )
        self._log("INFO", f"正在读取 {Path(path).name}")

        self.load_worker = LoadWorker(path, self.group_file, self)
        self.load_worker.finished_ok.connect(lambda ds, p=path: self._on_dataset_ready(ds, p, quiet))
        self.load_worker.failed.connect(self._on_dataset_failed)
        self.load_worker.start()

    def _on_dataset_ready(self, dataset, path: str, quiet: bool) -> None:
        self.dataset = dataset
        self.dataset_path = Path(path)
        self.file_badge.set_state("ok", "已加载")
        self.import_badge.set_state("ok", "已加载")
        self.dataset_summary.setText(dataset.summary())
        self.import_path_label.setText(f"{dataset.source}　·　{dataset.summary()}")
        self.output_dataset_label.setText(f"{dataset.source}\n{dataset.summary()}")
        self.top_badge.set_state("ok", "数据就绪")

        groups = dataset.group_names()
        for combo in (self.group_a_combo, self.group_b_combo):
            combo.clear()
            combo.addItems(groups)
        if len(groups) >= 2:
            self.group_a_combo.setCurrentIndex(0)
            self.group_b_combo.setCurrentIndex(1)
        self._log("SUCCESS", f"数据加载完成：{dataset.summary()}")
        for note in dataset.notes:
            self._log("INFO", note)
        self._fill_preview(dataset)
        self._apply_state()
        if not quiet:
            self._go(PAGE_ONE_CLICK)

    def _on_dataset_failed(self, message: str) -> None:
        self.dataset = None
        self.file_badge.set_state("err", "加载失败")
        self.import_badge.set_state("err", "加载失败")
        self.top_badge.set_state("err", "数据异常")
        self._log("ERROR", f"数据加载失败：{message}")
        self._apply_state()
        QMessageBox.critical(self, "数据加载失败", message)

    def _fill_preview(self, dataset) -> None:
        expr = dataset.expr
        head = expr.iloc[:12, :8]
        self.preview_table.clear()
        self.preview_table.setRowCount(head.shape[0])
        self.preview_table.setColumnCount(head.shape[1] + 1)
        headers = ["基因"] + [str(c) for c in head.columns]
        self.preview_table.setHorizontalHeaderLabels(headers)
        for r in range(head.shape[0]):
            item = QTableWidgetItem(str(head.index[r]))
            self.preview_table.setItem(r, 0, item)
            for c in range(head.shape[1]):
                value = head.iloc[r, c]
                cell = QTableWidgetItem(f"{value:.2f}" if isinstance(value, float) else str(value))
                cell.setTextAlignment(Qt.AlignRight | Qt.AlignVCenter)
                self.preview_table.setItem(r, c + 1, cell)
        self.preview_table.resizeColumnsToContents()
        self.preview_info.setText(
            f"共 {dataset.n_genes:,} 基因 × {dataset.n_samples:,} 样本　·　显示前 {head.shape[0]} 行 / {head.shape[1]} 列"
        )

    # ==================== 下载 ====================
    def _download_dialog(self) -> None:
        self._go(PAGE_IMPORT)
        self.url_edit.setFocus()

    def _download_from_edit(self) -> None:
        url = self.url_edit.text().strip()
        if not url:
            QMessageBox.information(self, "请输入链接", "请先粘贴数据文件的下载链接。")
            return
        if self.download_worker is not None and self.download_worker.isRunning():
            return
        dest = Path(self.output_edit.text().strip() or Path.home() / "OmicsOne_Data") / "downloads"
        self.url_go_btn.setEnabled(False)
        self.url_hint.setText("正在连接…")
        self.url_progress.setValue(0)

        self.download_worker = DownloadWorker(url, str(dest), self)
        self.download_worker.progress.connect(self._on_download_progress)
        self.download_worker.message.connect(lambda m: self.url_hint.setText(m))
        self.download_worker.finished_ok.connect(self._on_download_done)
        self.download_worker.failed.connect(self._on_download_failed)
        self.download_worker.start()

    def _on_download_progress(self, done: int, total: int) -> None:
        if total > 0:
            fraction = done / total
            self.url_progress.setValue(int(fraction * 1000))
            self.url_hint.setText(f"下载中 {done / 1e6:.1f} / {total / 1e6:.1f} MB（{fraction * 100:.0f}%）")
        else:
            self.url_hint.setText(f"下载中 已接收 {done / 1e6:.1f} MB")

    def _on_download_done(self, candidates) -> None:
        self.url_go_btn.setEnabled(True)
        self.url_progress.setValue(1000)
        paths = [Path(p) for p in candidates]
        if not paths:
            self._on_download_failed("未能识别下载内容")
            return
        self.url_hint.setText(f"下载完成，识别到 {len(paths)} 个候选文件")
        self._log("SUCCESS", f"下载完成：{paths[0].name}")
        if len(paths) == 1:
            self._load_file(str(paths[0]))
        else:
            choice, ok = QInputDialog.getItem(
                self, "选择数据文件", "压缩包里检测到多个数据文件：", [p.name for p in paths], 0, False
            )
            if ok and choice:
                self._load_file(str(paths[[p.name for p in paths].index(choice)]))

    def _on_download_failed(self, message: str) -> None:
        self.url_go_btn.setEnabled(True)
        self.url_hint.setText(f"下载失败：{message}")
        self._log("ERROR", f"下载失败：{message}")
        QMessageBox.critical(self, "下载失败", message)

    # ==================== 运行 ====================
    def _selected_steps(self) -> list[str]:
        steps = []
        for key, _label, step, _tip in ANALYSIS_ITEMS:
            if self.checkboxes[key].isChecked() and step not in steps:
                steps.append(step)
        return steps

    def _base_config(self) -> AnalysisConfig:
        output = self.output_edit.text().strip()
        if not output:
            output = str(Path.home() / "OmicsOne_Results")
            self.output_edit.setText(output)
            self.output_label.setText(output)
        name = self.project_name_edit.text().strip() if hasattr(self, "project_name_edit") else ""
        return AnalysisConfig(
            output_dir=Path(output),
            project_name=name if name and name != "未命名项目" else None,
            padj_cut=self.padj_spin.value(),
            lfc_cut=self.lfc_spin.value(),
            normalization=self.norm_combo.currentData(),
            gsea_permutations=self.perm_spin.value(),
            run_wgcna=self.checkboxes["wgcna"].isChecked(),
        )

    def _run_one_click(self) -> None:
        if self.dataset is None:
            return
        steps = self._selected_steps()
        if not steps:
            QMessageBox.information(self, "未选择分析项", "请至少勾选一项分析。")
            return
        config = self._base_config()

        ora_libs = []
        if self.checkboxes["kegg"].isChecked():
            ora_libs.append("KEGG_2021_Human")
        if self.checkboxes["go"].isChecked():
            ora_libs.append("GO_Biological_Process_2023")
        config.ora_libraries = ora_libs or ["KEGG_2021_Human"]
        config.gsea_library = "KEGG_2021_Human"

        if self.dataset.groups is None and any(s in steps for s in ("deg", "enrichment", "gsea")):
            QMessageBox.warning(
                self,
                "缺少分组信息",
                "该数据集没有识别到分组，差异分析、富集和 GSEA 需要分组。\n\n"
                "点「自动生成分组表」可按样本名推断分组并导出 CSV；"
                "也可以手工准备一张两列分组表（样本名, 组名），用「分组表（可选）」导入。",
            )
            return

        groups = self.dataset.group_names()
        if len(groups) >= 2:
            config.group_a = self.group_a_combo.currentText() or groups[0]
            config.group_b = self.group_b_combo.currentText() or groups[1]
            if config.group_a == config.group_b:
                QMessageBox.warning(self, "比较组无效", "比较组 A 与 B 不能相同。")
                return

        self._start_analysis(config, steps, label="一键全流程分析")

    def _start_analysis(self, config: AnalysisConfig, steps: list[str], *, label: str) -> None:
        if self.worker is not None and self.worker.isRunning():
            return
        self._log("INFO", f"开始{label}：{'、'.join(STEP_LABELS.get(s, s) for s in steps)}")
        self.top_badge.set_state("warn", "运行中")
        self.task_label.setText(f"当前任务：{label} · 准备中")
        self.percent_label.setText("0%")
        self.eta_label.setText("ETA --:--")
        self.progress_bar.setValue(0)
        self.result_label.setText(f"{label} 正在运行…")
        self._started_at = time.time()
        self._last_fraction = 0.0
        self._current_task = f"{label} · 准备中"
        self.pause_btn.setText("暂停")
        self._eta_timer.start()

        self.worker = AnalysisWorker(self.dataset, config, steps, self)
        self.worker.progress.connect(self._on_progress)
        self.worker.log.connect(self._on_worker_log)
        self.worker.finished_ok.connect(self._on_finished)
        self.worker.failed.connect(self._on_failed)
        self.worker.finished.connect(self._apply_state)
        self.worker.start()
        self._apply_state()

    def _on_progress(self, fraction: float, message: str) -> None:
        self._last_fraction = max(self._last_fraction, fraction)
        self.progress_bar.setValue(int(self._last_fraction * 1000))
        self.percent_label.setText(f"{self._last_fraction * 100:.0f}%")
        self._current_task = message
        self.task_label.setText(f"当前任务：{message}")

    def _tick_eta(self) -> None:
        elapsed = time.time() - self._started_at
        if self._last_fraction <= 0.01:
            return
        remaining = elapsed / self._last_fraction * (1 - self._last_fraction)
        minutes, seconds = divmod(int(max(remaining, 0)), 60)
        self.eta_label.setText(f"ETA {minutes:02d}:{seconds:02d}")

    def _on_worker_log(self, level: str, message: str) -> None:
        self._log(level, message)

    def _on_finished(self, result, report: str) -> None:
        self._eta_timer.stop()
        self.pause_btn.setText("暂停")
        self._paused_at = 0.0
        self.last_result = result
        self.last_report = report
        self.progress_bar.setValue(1000)
        self.percent_label.setText("100%")
        self.eta_label.setText("ETA 00:00")
        self.task_label.setText("当前任务：已完成")
        self.top_badge.set_state("ok", "已完成")
        self.result_label.setText(f"结果已写入：{result.project_dir}")
        self._log("SUCCESS", f"分析完成，结果目录 {result.project_dir}")
        self._log("INFO", f"HTML 报告：{Path(report).name}")
        if result.summary.get("elapsed_seconds"):
            self._log("INFO", f"总用时 {result.summary['elapsed_seconds']} 秒，产出图表 {len(result.figures)} 张")
        if result.errors:
            for err in result.errors:
                self._log("WARNING", f"未完成：{err}")
        self._apply_state()
        if result.errors and len(result.errors) == 1:
            self._log("INFO", "可在「运行日志」页查看完整堆栈")

    def _on_failed(self, message: str) -> None:
        self._eta_timer.stop()
        self.pause_btn.setText("暂停")
        self._paused_at = 0.0
        self.top_badge.set_state("err", "已终止" if message == "已终止" else "失败")
        self.task_label.setText("当前任务：已停止")
        self._log("ERROR", message if message != "已终止" else "分析已被用户终止")
        self._apply_state()

    def _toggle_pause(self) -> None:
        worker = self.worker
        if worker is None or not worker.isRunning():
            return
        if self.pause_btn.text() == "暂停":
            worker.pause()
            self._paused_at = time.time()
            self.pause_btn.setText("继续")
            self.top_badge.set_state("warn", "已暂停")
            self.task_label.setText(f"当前任务：{self._current_task} · 已暂停")
            self._log("WARNING", "已暂停：当前步骤跑完后停在检查点，点「继续」恢复")
        else:
            worker.resume()
            if self._paused_at:
                self._started_at += time.time() - self._paused_at  # 暂停时长不计入 ETA
                self._paused_at = 0.0
            self.pause_btn.setText("暂停")
            self.top_badge.set_state("warn", "运行中")
            self._log("INFO", "已继续")

    def _stop_run(self) -> None:
        if self.worker is not None and self.worker.isRunning():
            self.worker.cancel()
            self._log("WARNING", "正在终止分析…")
            self.stop_btn.setEnabled(False)

    # ==================== 打开结果 ====================
    def _choose_output_dir(self) -> None:
        folder = QFileDialog.getExistingDirectory(self, "选择结果输出文件夹", self.output_edit.text() or str(Path.home()))
        if folder:
            self.output_edit.setText(folder)
            self.output_label.setText(folder)
            self._log("INFO", f"输出目录：{folder}")

    def _open_report(self) -> None:
        if self.last_report:
            QDesktopServices.openUrl(QUrl.fromLocalFile(self.last_report))

    def _open_result_folder(self) -> None:
        if self.last_result is not None:
            open_folder(self.last_result.project_dir)

    def _export_log(self) -> None:
        text = self.log_panel.text() + "\n" + self.one_click_log.text()
        path, _ = QFileDialog.getSaveFileName(self, "导出日志", "OmicsOne_log.txt", "文本文件 (*.txt)")
        if path:
            Path(path).write_text(text, encoding="utf-8")
            self._log("INFO", f"日志已导出：{path}")

    def _show_help(self) -> None:
        QMessageBox.information(
            self,
            "使用说明",
            "1. 用「本地上传 / 链接下载 / 选择文件夹」导入表达矩阵（基因 × 样本，或样本 × 基因，自动识别）。\n"
            "2. 若样本名里没有组别前缀，点「自动生成分组表」按样本名推断分组并导出 CSV；也可以导入一张两列分组表：样本名, 组名。\n"
            "3. 勾选要跑的分析，设置输出文件夹，点击「一键开始分析」。\n"
            "4. 完成后点「打开报告」查看自包含 HTML。\n\n"
            "支持格式：csv / tsv / txt / tab / xlsx / xls / csv.gz；支持从链接下载 zip / tar.gz 后自动解压。\n"
            "全部计算在本机完成。",
        )

    # ==================== 日志 ====================
    def _log(self, level: str, message: str) -> None:
        stamp = time.strftime("%H:%M:%S")
        for panel in (getattr(self, "one_click_log", None), getattr(self, "log_panel", None)):
            if panel is not None:
                panel.append(level, message, stamp)

    def closeEvent(self, event) -> None:  # noqa: N802
        if self.worker is not None and self.worker.isRunning():
            answer = QMessageBox.question(self, "分析进行中", "分析仍在运行，确定退出吗？")
            if answer != QMessageBox.Yes:
                event.ignore()
                return
            self.worker.cancel()
            self.worker.wait(3000)
        event.accept()
