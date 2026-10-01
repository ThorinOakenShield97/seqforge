import typer

from seqforge.commands.input import InputSource, resolve_input
from seqforge.models.fastq_read import FastqRead
from seqforge.exceptions import InvalidFastaError
from seqforge.models.sequence import Sequence

def kmer(sequence: str, k:int | None = None, find: list[str] | None = None, counts: bool = False, frequencies: bool = False):
    """Display k-mers found in a biological sequence.

    Args:
        sequence: Literal sequence or path to a FASTA/FASTQ file.
        k: Length of each k-mer.
        counts: If True, display the counts of each k-mer.
     
     Raises:
          ValueError: If k is not a positive integer.
          FileNotFoundError: If the input file does not exist.
    """
    try:
        if counts and frequencies or counts and find or frequencies and find:
             raise ValueError('commands must be given one by one')
      
        results = resolve_input(sequence)
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
                   
    except FileNotFoundError as f:
            typer.echo(f"Error: {f}", err = True)
            raise typer.Exit(code=1)
    