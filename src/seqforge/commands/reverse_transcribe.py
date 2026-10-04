import typer

from seqforge.commands.input import resolve_input
from seqforge.models.sequence import Sequence
from seqforge.models.molecule_type import MoleculeType
from seqforge.exceptions import InvalidFastaError, InvalidFastqError


def reverse_transcribe(sequence:str) -> None:
    """
    Convert an RNA sequence into complementary DNA.

    Args:
        sequence: RNA sequence or path to a FASTA or FASTQ file.

    Raises:
        InvalidFastaError: If the FASTA input is invalid.
        InvalidFastqError: If the FASTQ input is invalid.
        FileNotFoundError: If the input file does not exist.
        ValueError: If the input sequence is not valid RNA.
    """

    try:
        results = resolve_input(sequence, molecule_type = MoleculeType.RNA)
    
        records = [Sequence(id = record.id, sequence = record.sequence, molecule_type = MoleculeType.RNA) for record in results.records]
        

        for record in records:
            reverse = record.reverse_transcribe()
            print(reverse.sequence)

    except InvalidFastaError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code=1)
    except InvalidFastqError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code=1)
    except FileNotFoundError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code=1)
    except ValueError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code=1)