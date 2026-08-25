from pathlib import Path

from seqforge.parsers.fastq import parse_fastq
from seqforge.exceptions import InvalidFastaError

import pytest



def test_parse_fastq_returns_read(tmp_path):
    fastq = tmp_path / "reads.fastq"
    fastq.write_text(
        "@read1\n"
        "ATGC\n"
        "+\n"
        "IIII\n"
    )

    reads = parse_fastq(fastq)

    assert len(reads) == 1
    assert reads[0].id == "read1"
    assert reads[0].sequence == "ATGC"
    assert reads[0].quality == "IIII"