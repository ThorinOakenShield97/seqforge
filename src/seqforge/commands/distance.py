import typer

from seqforge.models.sequence import Sequence

def distance(seq_1:str, seq_2:str, molecule_type:str = 'DNA') -> None:
    """
    Calculate the Hamming distance between two biological sequences.

    Args:
        seq_1: First sequence.
        seq_2: Second sequence.
        molecule_type: Molecule type of both input sequences. Defaults to DNA.

    Raises:
        ValueError: If the molecule type is invalid, the sequences have
        different lengths, or their molecule types are incompatible.
    """
    try:
        seq_1 = Sequence(id = 'seq 1', sequence = seq_1, molecule_type = molecule_type)
        seq_2 = Sequence(id = 'seq 2', sequence = seq_2, molecule_type = molecule_type)

        result = seq_1.distance(seq_2)

        print(f"distance: {result}")

    except ValueError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code=1)