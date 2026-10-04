import typer

from seqforge.models.sequence import Sequence
from seqforge.commands.input import InputSource, resolve_input
from seqforge.exceptions import InvalidFastaError, InvalidFastqError


def translate(sequence: str, frame: int | None = None, molecule_type: str = 'dna') -> None:
    """Translate a DNA or RNA sequence into a protein sequence.

    Args:
        sequence: Literal sequence or path to a FASTA or FASTQ file.
        frame: Reading frame to use, or None to search from the first start codon.
        molecule_type: Molecule type of the input sequence. Defaults to DNA.

    Raises:
        FileNotFoundError: If the input file does not exist.
        InvalidFastaError: If the FASTA input is invalid.
        InvalidFastqError: If the FASTQ input is invalid.
        ValueError: If the molecule type, frame, or sequence is invalid.
    """
    try:
        results = resolve_input(sequence, molecule_type)

        if results.source == InputSource.FASTQ:
            records = [Sequence(id=record.id, sequence=record.sequence, molecule_type=molecule_type) for record in results.records
    ]
        else:
            records = results.sequences

        for seq in records:
            if results.source == InputSource.FASTA:
                print(f">{seq.id}")
            elif results.source == InputSource.FASTQ:
                print(f"@{seq.id}")

            protein = seq.translate(frame = frame)
                                               
            if protein:
                print(f"Protein: {protein}")
            else:
                print("No start codon found.")
    
    except FileNotFoundError as f:
        typer.echo(f"Error: {f}", err=True)
        raise typer.Exit(code=1)
    except InvalidFastaError as i:
        typer.echo(f"Error: {i}", err=True)
        raise typer.Exit(code=1)
    except InvalidFastqError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code=1)
    except ValueError as v:
        typer.echo(f"Error: {v}", err=True)
        raise typer.Exit(code=1)
