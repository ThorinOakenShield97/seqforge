import typer

from seqforge.commands.input import InputSource, resolve_input
from seqforge.models.sequence import Sequence
from seqforge.exceptions import InvalidFastaError, InvalidFastqError


def motif(sequence:str, pattern:str, molecule_type:str = 'DNA') -> None:
    """
    Find all occurrences of a motif in a biological sequence.

    Args:
        sequence: Literal sequence or path to a FASTA or FASTQ file.
        pattern: Motif to search for.
        molecule_type: Molecule type of the input sequence. Defaults to DNA.

    Raises:
        InvalidFastaError: If the FASTA input is invalid.
        InvalidFastqError: If the FASTQ input is invalid.
        FileNotFoundError: If the input file does not exist.
        ValueError: If the molecule type, pattern, or sequence data is invalid.
    """

    try:
        results = resolve_input(sequence, molecule_type)

        if results.source == InputSource.FASTQ:
            records = [Sequence(id = record.id, sequence = record.sequence, molecule_type = molecule_type) for record in results.records]
    
        else:
            records = results.records

        for record in records:
            motifs = record.find_motif(pattern)
            print(f"positions: {motifs}")

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