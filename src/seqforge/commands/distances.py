import typer

from seqforge.commands.input import resolve_input, InputSource
from seqforge.models.sequence import Sequence
from seqforge.models.distances import pairwise_distances
from seqforge.exceptions import InvalidFastaError, InvalidFastqError

def distances(sequence:str, molecule_type:str = 'DNA') -> None:
    """
    Calculate pairwise Hamming distances for multiple sequences.

    Args:
        sequence: Literal sequence or path to a FASTA or FASTQ file.
        molecule_type: Molecule type shared by all input sequences. Defaults to DNA.

    Raises:
        InvalidFastaError: If the FASTA input is invalid.
        InvalidFastqError: If the FASTQ input is invalid.
        FileNotFoundError: If the input file does not exist.
        ValueError: If the molecule type or sequence data is invalid.
    """

    try:
        results = resolve_input(sequence, molecule_type = molecule_type)

        if results.source == InputSource.FASTQ:
            records = [Sequence(id = record.id, sequence = record.sequence, molecule_type = molecule_type) for record in results.records]

        else:
            records = results.records
    
        distances = pairwise_distances(records)

        for seq_1, seq_2, distance in distances:
            print(f"{seq_1}, {seq_2}: {distance}")

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