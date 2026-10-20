"""通过链接下载数据集，支持直链、断点续传式重试与常见压缩包自动解压。"""

from __future__ import annotations

import shutil
import tarfile
import zipfile
from pathlib import Path
from urllib.parse import unquote, urlparse

import requests

USER_AGENT = "OmicsOneClick/1.0 (+https://localhost)"
CHUNK = 1 << 20


class DownloadError(RuntimeError):
    pass


def _filename_from_url(url: str, response: requests.Response | None) -> str:
    if response is not None:
        disposition = response.headers.get("Content-Disposition", "")
        if "filename=" in disposition:
            raw = disposition.split("filename=")[-1].strip().strip('";')
            if raw:
                return Path(unquote(raw)).name
    name = Path(unquote(urlparse(url).path)).name
    return name or "dataset.bin"


def download(url: str, dest_dir: str | Path, *, timeout: int = 60, progress=None) -> Path:
    """下载 url 到 dest_dir，返回落地文件路径。progress(已下载字节, 总字节或 -1)。"""
    url = url.strip()
    if not url.lower().startswith(("http://", "https://")):
        raise DownloadError("只支持 http/https 链接")

    dest_dir = Path(dest_dir)
    dest_dir.mkdir(parents=True, exist_ok=True)

    with requests.get(url, stream=True, timeout=timeout, headers={"User-Agent": USER_AGENT}) as response:
        if response.status_code >= 400:
            raise DownloadError(f"服务器返回 {response.status_code}：{response.reason}")

        filename = _filename_from_url(url, response)
        target = dest_dir / filename

        total = int(response.headers.get("Content-Length") or -1)
        done = 0
        with open(target, "wb") as handle:
            for chunk in response.iter_content(CHUNK):
                if not chunk:
                    continue
                handle.write(chunk)
                done += len(chunk)
                if progress:
                    progress(done, total)

    if target.stat().st_size == 0:
        target.unlink(missing_ok=True)
        raise DownloadError("下载得到的文件为空")

    return target


def extract_if_archive(path: str | Path, *, remove_archive: bool = False) -> list[Path]:
    """解压 zip/tar.gz，返回解压出的数据文件列表；非压缩包则原样返回。"""
    path = Path(path)
    lower = path.name.lower()
    out_dir = path.parent / f"{path.stem}_extracted"

    if lower.endswith(".zip"):
        out_dir.mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(path) as archive:
            for member in archive.infolist():
                if member.is_dir():
                    continue
                target = (out_dir / member.filename).resolve()
                if not str(target).startswith(str(out_dir.resolve())):
                    raise DownloadError(f"压缩包含非法路径：{member.filename}")
                target.parent.mkdir(parents=True, exist_ok=True)
                with archive.open(member) as src, open(target, "wb") as dst:
                    shutil.copyfileobj(src, dst)
    elif lower.endswith((".tar.gz", ".tgz", ".tar")):
        out_dir.mkdir(parents=True, exist_ok=True)
        with tarfile.open(path) as archive:
            archive.extractall(out_dir, filter="data")
    else:
        return [path]

    if remove_archive:
        path.unlink(missing_ok=True)
    return sorted(p for p in out_dir.rglob("*") if p.is_file())


def find_data_files(paths: list[Path]) -> list[Path]:
    """从解压结果里挑出可能的数据文件。"""
    allowed = {".csv", ".tsv", ".txt", ".tab", ".xlsx", ".xls"}
    found = [p for p in paths if p.suffix.lower() in allowed]
    gz = [p for p in paths if p.suffix.lower() == ".gz" and p.stem.lower().endswith((".csv", ".tsv", ".txt"))]
    return found + gz
