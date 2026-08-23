from pathlib import Path

import typer

from seqforge.models.sequence import Sequence
from seqforge.parsers.fasta import parse_fasta
from seqforge.exceptions import InvalidFastaError

FASTA_EXTENSIONS = {".fasta", ".fa", ".fna"}

def stats(sequence: str) -> None:
    """Display sequence statistics."""
    if not sequence:
        typer.echo("Error: Cannot calculate statistics for an empty sequence.", err=True)
        raise typer.Exit(code=1)

    input_path = Path(sequence)

    if input_path.is_file():
        try:
            sequences = parse_fasta(input_path)
        except InvalidFastaError as e:
            typer.echo(f"Error: {e}", err=True)
            raise typer.Exit(code=1)
        show_id = True
    elif input_path.suffix.lower() in FASTA_EXTENSIONS:
        typer.echo(f"Error: Input file not found: {sequence}", err=True)
        raise typer.Exit(code=1)
    else:
        sequences = [Sequence(id="cli", sequence=sequence)]
        show_id = False



    try:
        for seq in sequences:
            if show_id:
                print(f">{seq.id}")

            counts = seq.base_counts()
            frequencies = seq.base_frequencies()

            print(f"Length: {seq.length()}")
            print(f"GC content: {seq.gc_content()}%")

            print("Base counts:")
            for base in sorted(counts):
                print(f"{base}: {counts[base]}")

            print("Base frequencies:")
            for base in sorted(frequencies):
                print(f"{base}: {frequencies[base]}%")
    except ValueError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code=1)