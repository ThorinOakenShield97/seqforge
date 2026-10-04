import typer

from seqforge.commands.version import version
from seqforge.commands.gc import gc
from seqforge.commands.translate import translate
from seqforge.commands.orf import orf
from seqforge.commands.stats import stats
from seqforge.commands.transcribe import transcribe
from seqforge.commands.kmer import kmer
from seqforge.commands.filter import filter
from seqforge.commands.align import align
from seqforge.commands.motif import motif
from seqforge.commands.reverse_transcribe import reverse_transcribe
from seqforge.commands.distance import distance
from seqforge.commands.distances import distances



app = typer.Typer(
    help="A modern toolkit for biological sequence analysis."
)


@app.callback()
def main() -> None:
    """SeqForge command line interface."""
    pass

app.command("align")(align)

app.command('distance')(distance)

app.command('distances')(distances)

app.command('filter')(filter)

app.command("gc")(gc)

app.command("kmer")(kmer)

app.command('motif')(motif)

app.command("orf")(orf)

app.command('reverse-transcribe')(reverse_transcribe)

app.command("stats")(stats)

app.command("transcribe")(transcribe)

app.command("translate")(translate)

app.command("version")(version)