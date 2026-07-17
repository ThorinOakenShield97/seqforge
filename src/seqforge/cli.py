import typer

app = typer.Typer(
    help="A modern toolkit for biological sequence analysis."
)


@app.callback()
def main() -> None:
    """SeqForge command line interface."""
    pass


@app.command()
def hello() -> None:
    """Test command."""

    print("Hello from SeqForge!")