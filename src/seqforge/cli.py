import typer

from seqforge.commands.version import version
from seqforge.commands.gc import gc
from seqforge.commands.translate import translate
from seqforge.commands.orf import orf
from seqforge.commands.stats import stats
from seqforge.commands.transcribe import transcribe
from seqforge.commands.kmer import kmer
from seqforge.commands.filter import filter




app = typer.Typer(
    help="A modern toolkit for biological sequence analysis."
)


@app.callback()
def main() -> None:
    """SeqForge command line interface."""
    pass

app.command("gc")(gc)

app.command("translate")(translate)

app.command("orf")(orf)

app.command("version")(version)

app.command("stats")(stats)

app.command("transcribe")(transcribe)

app.command("kmer")(kmer)

app.command('filter')(filter)