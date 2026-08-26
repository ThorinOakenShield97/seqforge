from pathlib import Path

from seqforge.parsers.fastq import parse_fastq
from seqforge.exceptions import InvalidFastqError
from seqforge.models.fastq_read import FastqRead
from seqforge.models.sequence import Sequence

import pytest

def test_fastq_read_returns_phred_scores():
    read = FastqRead(
        id="read1",
        sequence="ATGC",
        quality="IIII",
    )

    assert read.quality_scores() == [40, 40, 40, 40]

def test_fastq_read_quality_scores_match_sequence_length():
    read = FastqRead(
        id="read1",
        sequence="ATGC",
        quality="IIII",
    )

    assert len(read.quality_scores()) == len(read.sequence)

def test_fastq_read_mean_quality():
    read = FastqRead(
        id="read1",
        sequence="ATGC",
        quality="IIII",
    )

    assert read.mean_quality() == 40.0

def test_fastq_read_min_quality():
    read = FastqRead(
        id="read1",
        sequence="ATGC",
        quality="IHGF",
    )

    assert read.min_quality() == 37

def test_fastq_read_max_quality():
    read = FastqRead(
        id="read1",
        sequence="ATGC",
        quality="IHGF",
    )

    assert read.max_quality() == 40

def test_fastq_read_rejects_quality_length_mismatch():
    with pytest.raises(ValueError):
        FastqRead(
            id="read1",
            sequence="ATGC",
            quality="III",
        )

def test_fastq_read_single_base_quality():
    read = FastqRead(
        id="read1",
        sequence="A",
        quality="I",
    )

    assert read.quality_scores() == [40]
    assert read.mean_quality() == 40.0
    assert read.min_quality() == 40
    assert read.max_quality() == 40

def test_fastq_read_rejects_empty_sequence():
    with pytest.raises(ValueError):
        FastqRead(
            id="read1",
            sequence="",
            quality="",
        )

def test_fastq_read_rejects_invalid_quality_character():
    with pytest.raises(ValueError):
        FastqRead(
            id="read1",
            sequence="A",
            quality="\x20",
        )

def test_fastq_read_sequence_supports_kmers():
    read = FastqRead(
        id="read1",
        sequence="ATAT",
        quality="IIII",
    )

    assert Sequence(
        id=read.id,
        sequence=read.sequence,
    ).kmers(k=2) == [
        "AT",
        "TA",
        "AT",
    ]
