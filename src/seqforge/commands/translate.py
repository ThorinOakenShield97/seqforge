import typer

from seqforge.models.sequence import Sequence


def translate(sequence: str, frame: int | None = None) -> None:
    """Translate a DNA sequence into a protein sequence."""
    protein = Sequence(id="cli", sequence=sequence).translate(frame=frame)

    if protein:
        typer.echo(f"Protein: {protein}")
    else:
        typer.echo("No start codon found.")