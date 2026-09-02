from pathlib import Path
from dataclasses import dataclass
from seqforge.exceptions import InvalidFastqError
from seqforge.models.fastq_read import FastqRead


def parse_fastq(path: Path) -> list[FastqRead]:
    """Parse a FASTQ file and return its reads.

    Args:
        path: Path to the FASTQ file.

    Returns:
        A list of FastqRead objects.

    Raises:
        InvalidFastqError: If the FASTQ format is invalid.
    """
    records = []

    with open(path, 'r') as file:
        for header in file:
            header = header.strip()
            current_seq = file.readline().strip()
            separator = file.readline().strip()
            current_quality = file.readline().strip()

            if not header.startswith('@'):
                raise InvalidFastqError('Invalid Header')

            if not header[1:].strip():
                raise InvalidFastqError('Empty id')
            
            if separator != '+':
                raise InvalidFastqError('Invalid FASTQ separator')

            if len(current_seq) != len(current_quality):
                raise InvalidFastqError('Invalid FASTQ format')

            my_seq = FastqRead(id = header[1:], sequence = current_seq, quality = current_quality)
            records.append(my_seq)
            
    return records