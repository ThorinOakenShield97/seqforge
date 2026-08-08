from seqforge.models.sequence import Sequence

def test_seq_length():
    seq = Sequence(
        id="seq1",
        sequence="ATGC"
    )

    assert seq.length() == 4