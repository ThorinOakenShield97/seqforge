import typer

from seqforge.models.sequence import Sequence,MoleculeType
from seqforge.commands.input import InputSource, resolve_input
from seqforge.exceptions import InvalidFastaError
from seqforge.models.fastq_read import FastqRead


def gc(sequence: str, molecule_type: str = 'dna') -> None:
    """Display the GC content of a DNA sequence."""
    if not sequence:
        typer.echo("Error: Cannot calculate GC content of an empty sequence.", err=True)
        raise typer.Exit(code=1)

    try:
        results = resolve_input(sequence, molecule_type)

        for record in results.records:
            if isinstance(record, FastqRead):
                seq = Sequence(id=record.id, sequence=record.sequence, molecule_type = molecule_type)
                print(f"@{record.id}")
                print(f"GC content: {seq.gc_content()}%")

            elif results.source == InputSource.FASTA:
                print(f">{record.id}")
                print(f"GC content: {record.gc_content()}%")

            else:
                print(f"GC content: {record.gc_content()}%")

    except ValueError as v:
        typer.echo(f"Error: {v}", err=True)
        raise typer.Exit(code=1)
    except FileNotFoundError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code=1)
    except InvalidFastaError as i:
        typer.echo(f"Error: {i}", err=True)
        raise typer.Exit(code=1)
