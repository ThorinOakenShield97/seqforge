import typer

from seqforge.commands.version import app as version_app

app = typer.Typer(
    help="A modern toolkit for biological sequence analysis."
)


@app.callback()
def main() -> None:
    """SeqForge command line interface."""
    pass


app.add_typer(
    version_app,
    name="version",
)