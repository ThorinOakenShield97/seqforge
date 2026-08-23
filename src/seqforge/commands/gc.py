import typer

from seqforge.models.sequence import Sequence


def gc(sequence: str) -> None:
    """Display the GC content of a DNA sequence."""
    try:
        result = Sequence(id="cli", sequence=sequence).gc_content()
    except ValueError as exc:
        typer.echo(f"Error: {exc}", err=True)
        raise typer.Exit(code=1)

    print(f"GC content: {result}%")