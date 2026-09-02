import typer

from seqforge.models.sequence import Sequence,MoleculeType
from seqforge.commands.input import InputSource, resolve_input
from seqforge.exceptions import InvalidFastaError


def translate(sequence: str, frame: int | None = None, molecule_type: str = 'dna') -> None:
    """Translate a DNA sequence into a protein sequence."""
    try:
        results = resolve_input(sequence,molecule_type)

        for seq in results.sequences:
            if results.source == InputSource.FASTA:
                print(f">{seq.id}")

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
    except ValueError as v:
        typer.echo(f"Error: {v}", err=True)
        raise typer.Exit(code=1)
