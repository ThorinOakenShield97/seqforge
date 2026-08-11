import platform

import typer

from importlib.metadata import version as package_version


app = typer.Typer()


@app.callback(invoke_without_command=True)
def version() -> None:
    """Display SeqForge version information."""
    print("SeqForge")
    print(f"Version: {package_version('seqforge')}")
    print(f"Python: {platform.python_version()}")