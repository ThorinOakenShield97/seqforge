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