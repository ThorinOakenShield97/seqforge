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