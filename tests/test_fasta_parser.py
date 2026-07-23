from pathlib import Path

from seqforge.parsers.fasta import parse_fasta


def test_parse_single_sequence():
    records = parse_fasta(Path("tests/data/single_sequence.fasta"))
    # There is exactly one sequence.
    # The identifier is seq1.
    # The sequence is ATGC.
    assert len(records) == 1
    assert records[0].id == "seq1"
    assert records[0].sequence == "ATGC"

def test_parse_multiple_sequences():
    records = parse_fasta(Path("tests/data/multiple_sequence.fasta"))

    # There are two sequences.
    # The first and second identifiers are seq1 and seq2.
    # The sequences are ATGC and TTTT.
    assert len(records) == 2
    assert records[0].id == "seq1"
    assert records[0].sequence == "ATGC"
    assert records[1].id == 'seq2'
    assert records[1].sequence == 'TTTT'

def test_parse_multiline_sequence():
    records = parse_fasta(Path("tests/data/multiline_sequence.fasta"))

    # There is exactly one sequence.
    # The identifier is seq1.
    # The sequence spans multiple lines.    
    # The sequence is ATGCGGGGTTTT.
    assert len(records) == 1
    assert records[0].id == "seq1"
    assert records[0].sequence == "ATGCGGGGTTTT"