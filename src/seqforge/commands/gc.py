import typer
from seqforge.models.sequence import Sequence

app = typer.Typer()

@app.callback(invoke_without_command=True)
def gc(sequence: str) -> None:
    """Display the GC content of a DNA sequence."""
    try:
        result = Sequence(id="cli", sequence=sequence).gc_content()
    except ValueError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code = 1)
    
    print(f"GC content: {result}%")