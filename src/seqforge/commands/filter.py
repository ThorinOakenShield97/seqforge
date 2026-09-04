import typer

from seqforge.commands.input import InputSource, resolve_input
from seqforge.models.fastq_read import FastqRead
from seqforge.models.filters import filter_by_length, filter_by_motif, filter_by_quality

def filter(sequence: str, min_length: int | None = None, max_length: int | None = None, motif : str | None = None, min_quality: int | None = None) -> None:
    """Filter sequences or FASTQ reads using the selected criteria.

    Args:
        sequence: Literal sequence or path to a FASTA/FASTQ file.
        min_length: Minimum accepted sequence length, inclusive.
        max_length: Maximum accepted sequence length, inclusive.
        motif: Motif that must be present in the sequence.
        min_quality: Minimum accepted mean quality for FASTQ reads.

    Raises:
        ValueError: If no filters are specified or a filter is incompatible
            with the input type.
    """

    if min_length is None and max_length is None and motif is None and min_quality is None:
        raise ValueError('No filters specified')

    results = resolve_input(sequence)

    filtered = results.records

    if min_length is not None or max_length is not None:
        filtered = filter_by_length(filtered, min_length = min_length, max_length = max_length)

    if motif is not None:
        filtered = filter_by_motif(filtered, motif = motif)

    if min_quality is not None:
        if results.source == InputSource.FASTQ:
            filtered = filter_by_quality(filtered, min_quality=min_quality)

        else:
            raise ValueError("Cannot apply quality filtering to non-FASTQ input")


    for record in filtered:
        if results.source == InputSource.FASTA:
            print(f">{record.id}")
            print(record.sequence)

        elif isinstance(record, FastqRead):
            print(f"@{record.id}")
            print(record.sequence)
            print("+")
            print(record.quality)

        else:
            print(record.sequence)


