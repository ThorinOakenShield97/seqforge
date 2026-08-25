from pathlib import Path

from seqforge.parsers.fastq import parse_fastq
from seqforge.exceptions import InvalidFastqError
from seqforge.models.fastq_read import FastqRead

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

def test_parse_fastq_rejects_quality_length_mismatch(tmp_path):
    fastq = tmp_path / "invalid.fastq"
    fastq.write_text(
        "@read1\n"
        "ATGC\n"
        "+\n"
        "III\n"
    )

    with pytest.raises(InvalidFastqError):
        parse_fastq(fastq)

def test_parse_fastq_rejects_invalid_separator(tmp_path):
    fastq = tmp_path / "invalid.fastq"
    fastq.write_text(
        "@read1\n"
        "ATGC\n"
        "-\n"
        "IIII\n"
    )

    with pytest.raises(InvalidFastqError):
        parse_fastq(fastq)

def test_parse_fastq_returns_multiple_reads(tmp_path):
    fastq = tmp_path / "reads.fastq"
    fastq.write_text(
        "@read1\n"
        "ATGC\n"
        "+\n"
        "IIII\n"
        "@read2\n"
        "GGTA\n"
        "+\n"
        "HHHH\n"
    )

    reads = parse_fastq(fastq)

    assert len(reads) == 2

    assert reads[0].id == "read1"
    assert reads[0].sequence == "ATGC"
    assert reads[0].quality == "IIII"

    assert reads[1].id == "read2"
    assert reads[1].sequence == "GGTA"
    assert reads[1].quality == "HHHH"

def test_parse_fastq_rejects_incomplete_record(tmp_path):
    fastq = tmp_path / "invalid.fastq"
    fastq.write_text(
        "@read1\n"
        "ATGC\n"
        "+\n"
    )

    with pytest.raises(InvalidFastqError):
        parse_fastq(fastq)

def test_parse_fastq_rejects_invalid_header(tmp_path):
    fastq = tmp_path / "invalid.fastq"
    fastq.write_text(
        "read1\n"
        "ATGC\n"
        "+\n"
        "IIII\n"
    )

    with pytest.raises(InvalidFastqError):
        parse_fastq(fastq)

def test_parse_fastq_rejects_empty_header_id(tmp_path):
    fastq = tmp_path / "invalid.fastq"
    fastq.write_text(
        "@\n"
        "ATGC\n"
        "+\n"
        "IIII\n"
    )

    with pytest.raises(InvalidFastqError):
        parse_fastq(fastq)

def test_parse_fastq_preserves_full_header_id(tmp_path):
    fastq = tmp_path / "reads.fastq"
    fastq.write_text(
        "@read1 description of the read\n"
        "ATGC\n"
        "+\n"
        "IIII\n"
    )

    reads = parse_fastq(fastq)

    assert reads[0].id == "read1 description of the read"

def test_parse_fastq_preserves_quality_for_multiple_reads(tmp_path):
    fastq = tmp_path / "reads.fastq"
    fastq.write_text(
        "@read1\n"
        "ATGC\n"
        "+\n"
        "IIII\n"
        "@read2\n"
        "GGTA\n"
        "+\n"
        "HHHH\n"
    )

    reads = parse_fastq(fastq)

    assert reads[0].mean_quality() == 40.0
    assert reads[1].mean_quality() == 39.0

