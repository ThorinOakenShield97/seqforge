import typer

from seqforge.commands.version import app as version_app
from seqforge.commands.gc import app as gc_app
from seqforge.commands.translate import translate


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

app.add_typer(
    gc_app,
    name="gc",
)

app.command("translate")(translate)

"""
app.add_typer(
    translate_app,
    name="translate",
)"""