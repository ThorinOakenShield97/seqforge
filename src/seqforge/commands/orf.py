import typer

from seqforge.models.sequence import Sequence


app = typer.Typer()


@app.callback(invoke_without_command=True)
def orf(sequence: str, strand: str = "forward", frame: int | None = None) -> None:
    """Find open reading frames in a DNA sequence."""

    try:
        orfs = Sequence(id = 'cli', sequence = sequence).find_orfs(strand = strand, frame = frame)

    except ValueError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code = 1)
        
    for elem in orfs:
        print(elem)
    