# -*- mode: python ; coding: utf-8 -*-
"""PyInstaller 打包配置：把 GMT 基因集库与图标一并打进 exe。"""

from pathlib import Path

root = Path(SPECPATH)

datas = [
    (str(root / "app" / "assets" / "gmt"), "assets/gmt"),
    (str(root / "app" / "assets" / "logo.png"), "assets"),
]

hiddenimports = [
    "openpyxl",
    "xlrd",
    "scipy.special",
    "scipy.stats",
    "scipy.cluster.hierarchy",
    "scipy.spatial.distance",
    "matplotlib.backends.backend_agg",
    "pandas.io.formats.excel",
]

a = Analysis(
    [str(root / "app" / "main.py")],
    pathex=[str(root)],
    binaries=[],
    datas=datas,
    hiddenimports=hiddenimports,
    hookspath=[],
    runtime_hooks=[],
    excludes=[
        "PySide6.Qt3DCore", "PySide6.Qt3DRender", "PySide6.QtBluetooth", "PySide6.QtCharts",
        "PySide6.QtDataVisualization", "PySide6.QtDesigner", "PySide6.QtMultimedia",
        "PySide6.QtMultimediaWidgets", "PySide6.QtOpenGL", "PySide6.QtPdf", "PySide6.QtQml",
        "PySide6.QtQuick", "PySide6.QtQuick3D", "PySide6.QtRemoteObjects", "PySide6.QtSensors",
        "PySide6.QtSerialPort", "PySide6.QtSpatialAudio", "PySide6.QtSql", "PySide6.QtTest",
        "PySide6.QtWebChannel", "PySide6.QtWebEngineCore", "PySide6.QtWebEngineWidgets",
        "PySide6.QtWebSockets", "PySide6.QtWebView", "PySide6.QtScxml", "PySide6.QtStateMachine",
        "tkinter", "PyQt5", "PyQt6", "IPython", "jedi", "pytest",
    ],
    noarchive=False,
)

pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name="OmicsOne",
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=False,
    console=False,
    disable_windowed_traceback=False,
    icon=str(root / "app" / "assets" / "logo.ico"),
)
