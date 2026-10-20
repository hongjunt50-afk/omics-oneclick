"""后台工作线程：分析任务与下载任务，均不阻塞界面。"""

from __future__ import annotations

import traceback
from pathlib import Path

from PySide6.QtCore import QThread, Signal

from ..core.data_loader import load_dataset
from ..core.downloader import DownloadError, download, extract_if_archive, find_data_files
from ..core.pipeline import AnalysisConfig, Cancelled, Runner
from ..core.report import build_report


class LoadWorker(QThread):
    """读取并解析数据集。"""

    finished_ok = Signal(object)   # Dataset
    failed = Signal(str)

    def __init__(self, path: str, group_file: str | None = None, parent=None):
        super().__init__(parent)
        self.path = path
        self.group_file = group_file

    def run(self) -> None:
        try:
            dataset = load_dataset(self.path, group_file=self.group_file)
            self.finished_ok.emit(dataset)
        except Exception as exc:
            self.failed.emit(f"{exc}")


class DownloadWorker(QThread):
    """下载链接 → 落地文件 → 必要时解压 → 返回候选数据文件。"""

    progress = Signal(int, int)      # 已下载, 总大小
    message = Signal(str)
    finished_ok = Signal(list)       # list[Path]
    failed = Signal(str)

    def __init__(self, url: str, dest: str, parent=None):
        super().__init__(parent)
        self.url = url
        self.dest = dest

    def run(self) -> None:
        try:
            self.message.emit("正在下载…")
            path = download(self.url, self.dest, progress=self.progress.emit)
            self.message.emit(f"下载完成：{path.name}")
            extracted = extract_if_archive(path)
            candidates = find_data_files(extracted) if len(extracted) > 1 else [path]
            if not candidates:
                candidates = find_data_files([path])
            if not candidates:
                raise DownloadError("下载的文件不是可识别的数据表（支持 csv/tsv/txt/xlsx 或包含它们的压缩包）")
            self.finished_ok.emit(candidates)
        except DownloadError as exc:
            self.failed.emit(str(exc))
        except Exception as exc:
            self.failed.emit(f"下载失败：{exc}")


class AnalysisWorker(QThread):
    """跑流水线 + 出 HTML 报告。"""

    progress = Signal(float, str)
    log = Signal(str, str)
    finished_ok = Signal(object, str)   # RunResult, report_path
    failed = Signal(str)

    def __init__(self, dataset, config: AnalysisConfig, steps: list[str] | None = None, parent=None):
        super().__init__(parent)
        self.dataset = dataset
        self.config = config
        self.steps = steps
        self._runner: Runner | None = None

    def cancel(self) -> None:
        if self._runner is not None:
            self._runner.cancel()
        self.requestInterruption()

    def pause(self) -> None:
        if self._runner is not None:
            self._runner.pause()

    def resume(self) -> None:
        if self._runner is not None:
            self._runner.resume()

    def run(self) -> None:
        try:
            self._runner = Runner(
                self.dataset,
                self.config,
                progress=lambda f, m: self.progress.emit(float(f), str(m)),
                log=lambda level, msg: self.log.emit(str(level), str(msg)),
            )
            result = self._runner.run_all(self.steps)
            self.progress.emit(0.95, "生成 HTML 报告")
            report = build_report(result, self.dataset, self.config)
            self.finished_ok.emit(result, str(report))
        except Cancelled:
            self.failed.emit("已终止")
        except Exception as exc:
            self.log.emit("ERROR", traceback.format_exc(limit=3))
            self.failed.emit(str(exc))
        finally:
            self._runner = None


def open_folder(path: str | Path) -> None:
    """在资源管理器里打开目录。"""
    import os
    import subprocess
    import sys

    target = str(path)
    if not Path(target).exists():
        Path(target).mkdir(parents=True, exist_ok=True)
    if sys.platform.startswith("win"):
        os.startfile(target)  # noqa: S606
    elif sys.platform == "darwin":
        subprocess.Popen(["open", target])
    else:
        subprocess.Popen(["xdg-open", target])
