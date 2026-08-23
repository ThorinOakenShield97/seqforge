
import typer

from seqforge.commands.input import InputSource, resolve_input
from seqforge.exceptions import InvalidFastaError


def stats(sequence: str) -> None:
    """Display sequence statistics."""
    try:
        resolved = resolve_input(sequence)
        if not sequence:
            typer.echo("Error: Cannot calculate statistics for an empty sequence.",err=True)
            raise typer.Exit(code=1)
        
        for seq in resolved.sequences:
            if resolved.source == InputSource.FASTA:
                print(f">{seq.id}")

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
    except FileNotFoundError as f:
        typer.echo(f"Error: {f}", err=True)
        raise typer.Exit(code=1)
    except InvalidFastaError as i:
        typer.echo(f"Error: {i}", err=True)
        raise typer.Exit(code=1)