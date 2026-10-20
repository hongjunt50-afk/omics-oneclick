"""GMT 基因集库的读取、归类与内置库清单。"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

# 内置库：文件名（不含 .gmt）-> 中文显示名
BUNDLED: dict[str, str] = {
    "KEGG_2021_Human": "KEGG 2021 人类通路",
    "GO_Biological_Process_2023": "GO 生物过程 (2023)",
    "GO_Molecular_Function_2023": "GO 分子功能 (2023)",
    "GO_Cellular_Component_2023": "GO 细胞组分 (2023)",
    "Reactome_2022": "Reactome 2022 通路",
    "MSigDB_Hallmark_2020": "MSigDB Hallmark 2020",
}


@dataclass
class GeneSetLibrary:
    name: str
    label: str
    sets: dict[str, list[str]]

    @property
    def size(self) -> int:
        return len(self.sets)


def asset_dir() -> Path:
    """打包后从 _MEIPASS 读，开发时从 app/assets 读。"""
    import sys

    base = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent.parent))
    candidates = [base / "assets" / "gmt", base / "gmt", Path(__file__).resolve().parent.parent / "assets" / "gmt"]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return candidates[0]


def read_gmt(path: str | Path, *, label: str | None = None) -> GeneSetLibrary:
    path = Path(path)
    sets: dict[str, list[str]] = {}
    with open(path, "r", encoding="utf-8", errors="replace") as handle:
        for line in handle:
            line = line.rstrip("\n").rstrip("\r")
            if not line.strip():
                continue
            fields = line.split("\t")
            if len(fields) < 3:
                continue
            name = fields[0].strip()
            genes = [g.strip().upper() for g in fields[2:] if g.strip()]
            if name and genes:
                sets[name] = genes
    if not sets:
        raise ValueError(f"GMT 文件里没有可用的基因集：{path.name}")
    return GeneSetLibrary(name=path.stem, label=label or path.stem, sets=sets)


def available_libraries(extra_dir: str | Path | None = None) -> dict[str, tuple[Path, str]]:
    """返回 库名 -> (路径, 显示名)，内置库 + 用户目录里的自定义 gmt。"""
    found: dict[str, tuple[Path, str]] = {}
    bundled_dir = asset_dir()
    for name, label in BUNDLED.items():
        path = bundled_dir / f"{name}.gmt"
        if path.exists():
            found[name] = (path, label)

    if extra_dir:
        extra = Path(extra_dir)
        if extra.exists():
            for path in sorted(extra.glob("*.gmt")):
                found.setdefault(path.stem, (path, path.stem))
    return found


def load_library(name: str, extra_dir: str | Path | None = None) -> GeneSetLibrary:
    libs = available_libraries(extra_dir)
    if name not in libs:
        raise KeyError(f"未找到基因集库 {name}")
    path, label = libs[name]
    return read_gmt(path, label=label)
