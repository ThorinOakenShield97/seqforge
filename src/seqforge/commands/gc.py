import typer

from seqforge.models.sequence import Sequence
from seqforge.commands.input import InputSource, resolve_input
from seqforge.exceptions import InvalidFastaError


def gc(sequence: str) -> None:
    """Display the GC content of a DNA sequence."""
    try:
        results = resolve_input(sequence)

        for seq in results.sequences:
            if results.source == InputSource.FASTA:
                print(f">{seq.id}")

            print(f"GC content: {seq.gc_content()}%")

    except ValueError as v:
        typer.echo(f"Error: {v}", err=True)
        raise typer.Exit(code=1)
    except FileNotFoundError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code=1)
    except InvalidFastaError as i:
        typer.echo(f"Error: {i}", err=True)
        raise typer.Exit(code=1)
