# -*- mode: python ; coding: utf-8 -*-
from PyInstaller.utils.hooks import collect_submodules

hiddenimports = ['matplotlib.backends.backend_tkagg']
hiddenimports += collect_submodules('pandas')


a = Analysis(
    ['C:\\Users\\eduso\\Documents\\GIT\\Python\\Rectifier_C_Project\\src\\rectifier_abacos_gui.py'],
    pathex=[],
    binaries=[],
    datas=[('C:\\Users\\eduso\\Documents\\GIT\\Python\\Rectifier_C_Project\\src\\rectifier_abacos_core.py', '.')],
    hiddenimports=hiddenimports,
    hookspath=[],
    hooksconfig={},
    runtime_hooks=[],
    excludes=['matplotlib.tests', 'matplotlib.backends._macosx', 'matplotlib.backends.backend_qt5', 'matplotlib.backends.backend_qt4', 'matplotlib.backends.backend_qtagg', 'matplotlib.backends.qt_compat', 'PyQt5', 'PySide2', 'PySide6', 'tkinter.test', 'numpy.tests', 'pandas.tests'],
    noarchive=False,
    optimize=0,
)
pyz = PYZ(a.pure)

exe = EXE(
    pyz,
    a.scripts,
    a.binaries,
    a.datas,
    [],
    name='RectifierAbacos',
    debug=False,
    bootloader_ignore_signals=False,
    strip=False,
    upx=True,
    upx_exclude=[],
    runtime_tmpdir=None,
    console=False,
    disable_windowed_traceback=False,
    argv_emulation=False,
    target_arch=None,
    codesign_identity=None,
    entitlements_file=None,
    icon=['C:\\Users\\eduso\\Documents\\GIT\\Python\\Rectifier_C_Project\\assets\\icon.ico'],
)
