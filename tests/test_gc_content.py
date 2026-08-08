from seqforge.models.sequence import Sequence

def test_gc_content():
    seq = Sequence(
        id="seq1",
        sequence="ATGC"
    )
    assert seq.gc_content() == 50.0