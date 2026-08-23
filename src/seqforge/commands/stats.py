import typer

from seqforge.models.sequence import Sequence

def stats(sequence: str) -> None:
    """Display sequence statistics."""
    if not sequence:
        typer.echo("Error: Cannot calculate statistics for an empty sequence.", err=True)
        raise typer.Exit(code=1)
    
    try:
        seq = Sequence(id="cli", sequence=sequence)
        counts = seq.base_counts()
        frequencies = seq.base_frequencies()

        print(f"Length: {seq.length()}")
        print(f"GC content: {seq.gc_content()}%")

        print("Base counts:")
        for base in sorted(counts):
            print(f"{base}: {counts[base]}")

        print("Base frequencies:")
        for base in sorted(frequencies):
            print(f"{base}: {frequencies[base]}%")

    except ValueError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code=1)
    