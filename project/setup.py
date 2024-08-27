from cx_Freeze import setup, Executable

excludes = [
    # List of packages/modules to exclude
    'tkinter',
    'PySide6',
]

build_options = {
    'excludes': excludes,
}

setup(
    name="FOME",
    version="1.0",
    description="Project",
    executables=[Executable("main.py")],
    options={'build_exe': build_options},
)