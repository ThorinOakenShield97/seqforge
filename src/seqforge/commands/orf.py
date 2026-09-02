from seqforge.commands.input import InputSource, resolve_input
from seqforge.exceptions import InvalidFastaError
from seqforge.models.sequence import Sequence, MoleculeType
import typer

def orf(sequence: str, strand: str = "forward", frame: int | None = None, molecule_type: str = 'dna') -> None:
    """Find open reading frames in a DNA or RNA sequence.

    Args:
        sequence: Literal sequence or path to a FASTA file.
        strand: Strand to search: "forward", "reverse", or "both".
        frame: Reading frame to search, or None to search all frames.
        molecule_type: Molecule type of the input sequence. Defaults to DNA.

    Raises:
        FileNotFoundError: If the input file does not exist.
        InvalidFastaError: If the FASTA file is invalid.
        ValueError: If the molecule type, strand, frame, or sequence is invalid."""

    try:
        results = resolve_input(sequence, molecule_type)

        for seq in results.sequences:
            if results.source == InputSource.FASTA:
                print(f">{seq.id}")

            
            if strand == 'both':
                forward_orfs = seq.find_orfs(strand = "forward", frame = frame)
                reverse_orfs = seq.find_orfs(strand = 'reverse', frame = frame)

                print("Forward:")
                if forward_orfs:
                    for forward in forward_orfs:
                        print(forward)
                else:
                    print("No ORFs found")

                print("Reverse:")
                if reverse_orfs:
                    for reverse in reverse_orfs:
                        print(reverse)
                else:
                    print("No ORFs found.")

            else:
                orfs = seq.find_orfs(strand = strand, frame = frame)

                for elem in orfs:
                    print(elem)

    except ValueError as v:
        typer.echo(f"Error: {v}", err=True)
        raise typer.Exit(code = 1)
    except InvalidFastaError as i:
        typer.echo(f"Error: {i}", err=True)
        raise typer.Exit(code=1)
    except FileNotFoundError as f:
        typer.echo(f"Error: {f}", err=True)
        raise typer.Exit(code=1)
        
    
    