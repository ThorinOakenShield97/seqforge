import typer

from seqforge.commands.input import InputSource, resolve_input
from seqforge.models.fastq_read import FastqRead
from seqforge.exceptions import InvalidFastaError
from seqforge.models.sequence import Sequence

def kmer(sequence: str, k:int | None = None, counts:bool = False):
    try:
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

            k_mers = seq.kmers(k)

            if counts:
                freqs = seq.kmer_counts(k)
                for k_mer in freqs:
                     print(f"{k_mer}: {freqs[k_mer]}")         
            else:
                for k_mer in k_mers:
                    print(k_mer)
                     
    except FileNotFoundError as f:
            typer.echo(f"Error: {f}", err=True)
            raise typer.Exit(code=1)
    