import typer

from seqforge.commands.input import InputSource, resolve_input
from seqforge.models.fastq_read import FastqRead
from seqforge.models.sequence import Sequence
from seqforge.exceptions import InvalidFastaError, InvalidFastqError

def kmer(sequence: str, k:int | None = None, find: list[str] | None = None, counts: bool = False, frequencies: bool = False, molecule_type: str = "DNA") -> None:
    """
    Display k-mers found in a biological sequence.

    Args:
        sequence: Literal sequence or path to a FASTA or FASTQ file.
        k: Length of each k-mer.
        find: K-mers to search for in the sequence.
        counts: If True, display the count of each k-mer.
        frequencies: If True, display the relative frequency of each k-mer.
        molecule_type: Molecule type of the input sequence. Defaults to DNA.

    Raises:
        InvalidFastaError: If the FASTA input is invalid.
        InvalidFastqError: If the FASTQ input is invalid.
        FileNotFoundError: If the input file does not exist.
        ValueError: If the k-mer parameters or molecule type are invalid,
        or if incompatible modes are selected.
    """
    
    try:
        if counts and frequencies or counts and find or frequencies and find:
             raise ValueError('commands must be given one by one')
      
        results = resolve_input(sequence, molecule_type)
        for record in results.records:
          if isinstance(record, FastqRead):
               seq = Sequence(id = record.id, sequence = record.sequence)
               print(f"@{record.id}")

          elif results.source == InputSource.FASTA:
               print(f">{record.id}")
               seq = record
          else:
               seq = record

          if counts:
               count = seq.kmer_counts(k)
               for k_mer in count:
                    print(f"{k_mer}: {count[k_mer]}")

          elif frequencies:
               freqs = seq.kmer_frequencies(k)
               for k_mer in freqs:
                    print(f"{k_mer}: {freqs[k_mer]}")

          elif find:
               found = seq.find_kmers(find)
               for k_mer in found:
                    print(f"{k_mer}: {found[k_mer]}")

          else:
               k_mers = seq.kmers(k)
               for k_mer in k_mers:
                    print(k_mer)
                   
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
    