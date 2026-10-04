import typer

from seqforge.commands.input import InputSource, resolve_input
from seqforge.models.sequence import Sequence
from seqforge.models.alignment import multiple_alignment
from seqforge.exceptions import InvalidFastqError, InvalidFastaError

def align(sequence:str, molecule_type: str = "DNA") -> None:
    """
    Align multiple biological sequences.

    Args:
        sequence: Literal sequence or path to a FASTA or FASTQ file.
        molecule_type: Molecule type shared by all input sequences.
            Defaults to DNA.

    Raises:
        InvalidFastaError: If the FASTA input is invalid.
        InvalidFastqError: If the FASTQ input is invalid.
        FileNotFoundError: If the input file does not exist.
        ValueError: If the molecule type or sequence data is invalid.
    """

    try:
        results = resolve_input(sequence, molecule_type)

        if results.source == InputSource.FASTQ:
            records = [Sequence(id = record.id, sequence = record.sequence, molecule_type = molecule_type) for record in results.records]

        else:
            records = results.records
    
        alignments = multiple_alignment(records)

        for alignment in alignments:
            print(alignment)

    except InvalidFastaError as e:
            typer.echo(f"Error: {e}", err=True)
            raise typer.Exit(code=1)
    except FileNotFoundError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code=1)
    except ValueError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code=1)
    except InvalidFastqError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code=1)