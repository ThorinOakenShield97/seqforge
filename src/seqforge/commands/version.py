import platform

import typer

app = typer.Typer()


@app.callback(invoke_without_command=True)
def version() -> None:
    """Display SeqForge version information."""

    print("SeqForge")
    print("Version: 0.1.0")
    print(f"Python: {platform.python_version()}")