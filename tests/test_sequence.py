import pytest

from seqforge.models.sequence import Sequence


def test_seq_length():
    seq = Sequence(
        id="seq1",
        sequence="ATGC"
    )

    assert seq.length() == 4

def test_gc_content_empty_seq():
    seq = Sequence(
        id="seq1",
        sequence=""
    )

    with pytest.raises(ValueError):
        seq.gc_content()

def test_reverse_complement():
    seq = Sequence(
        id="seq1",
        sequence="ATGC"
    )

    assert seq.reverse_complement() == "GCAT"


def test_reverse_complement_lowercase():
    seq = Sequence(
        id="seq1",
        sequence="atgc"
    )

    assert seq.reverse_complement() == "GCAT"

def test_reverse_complement_with_ambiguous_base():
    seq = Sequence(
        id="seq1",
        sequence="ATGN"
    )

    assert seq.reverse_complement() == "NCAT"

def test_reverse_complement_iupac():
    seq = Sequence(
        id="seq1",
        sequence="ACGTRYSWKMBDHVN"
    )

    assert seq.reverse_complement() == "NBDHVKMWSRYACGT"

def test_reverse_complement_invalid_base():
    seq = Sequence(
        id="seq1",
        sequence="ATGX"
    )

    with pytest.raises(ValueError):
        seq.reverse_complement()

def test_reverse_complement_empty_sequence():
    seq = Sequence(
        id="seq1",
        sequence=""
    )

    assert seq.reverse_complement() == ""