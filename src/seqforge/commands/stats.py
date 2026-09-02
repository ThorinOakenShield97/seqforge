
import typer

from seqforge.commands.input import InputSource, resolve_input
from seqforge.models.fastq_read import FastqRead
from seqforge.exceptions import InvalidFastaError
from seqforge.models.sequence import Sequence, MoleculeType


def stats(sequence: str, molecule_type: str = 'dna' ) -> None:
    """Display sequence statistics."""
    try:
        resolved = resolve_input(sequence, molecule_type)
        if not sequence:
            typer.echo("Error: Cannot calculate statistics for an empty sequence.",err=True)
            raise typer.Exit(code=1)

        reads = 0
        total_length = 0
        mean = 0
        gc = 0
        results = []
        for record in resolved.records:
            length = 0
            if isinstance(record, FastqRead):
                seq = Sequence(id=record.id, sequence=record.sequence)
                reads += 1
                total_length += seq.length()
                length += seq.length()
                results.append(length)
                mean += record.mean_quality()
                gc += seq.gc_content()
                
                print(f"@{record.id}")
                print(f"Length: {seq.length()}")
                print(f"GC content: {seq.gc_content()}%")
                print(f"Mean quality: {record.mean_quality()}")
                print(f"Min quality: {record.min_quality()}")
                print(f"Max quality: {record.max_quality()}")

            else:
                seq = record

                if resolved.source == InputSource.FASTA:
                    print(f">{seq.id}")

                if seq.molecule_type == MoleculeType.DNA or seq.molecule_type == MoleculeType.RNA:
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

                else:
                    counts = seq.amino_acid_counts()
                    frequencies = seq.amino_acid_frequencies()

                    print(f"Length: {seq.length()}")

                    print("Amino acid counts:")
                    for aa in sorted(counts):
                        print(f"{aa}: {counts[aa]}")

                    print("Amino acid frequencies:")
                    for aa in sorted(frequencies):
                        print(f"{aa}: {frequencies[aa]}%")


        if resolved.source == InputSource.FASTQ:
            print(f"Reads: {reads}")
            print(f"Mean read length: {total_length/reads}")
            print(f"Overall mean quality: {mean/reads}")
            print(f"Mean GC content: {gc/reads}%")
            print(f"Min read length: {min(results)}")
            print(f"Max read length: {max(results)}")

    except ValueError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code=1)
    except FileNotFoundError as f:
        typer.echo(f"Error: {f}", err=True)
        raise typer.Exit(code=1)
    except InvalidFastaError as i:
        typer.echo(f"Error: {i}", err=True)
        raise typer.Exit(code=1)