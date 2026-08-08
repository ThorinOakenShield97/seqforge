import pytest

from seqforge.models.sequence import Sequence, expand_iupac


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

def test_transcribe():
    seq = Sequence(
        id="seq1",
        sequence="ATGC"
    )

    assert seq.transcribe() == "AUGC"


def test_transcribe_lowercase():
    seq = Sequence(
        id="seq1",
        sequence="atgc"
    )

    assert seq.transcribe() == "AUGC"

def test_transcribe_iupac():
    seq = Sequence(
        id="seq1",
        sequence="ATGN"
    )

    assert seq.transcribe() == "AUGN"

def test_translate():
    seq = Sequence(
        id="seq1",
        sequence="ATGGCC"
    )

    assert seq.translate() == "MA"

def test_translate_stop():
    seq = Sequence(
        id="seq1",
        sequence="ATGGCCTAA"
    )

    assert seq.translate() == "MA"

def test_translate_without_start_codon():
    seq = Sequence(
        id="seq1",
        sequence="CCCGGGAAC"
    )

    assert seq.translate() == ""

def test_translate_incomplete_codon():
    seq = Sequence(
        id="seq1",
        sequence="ATGGCCA"
    )

    assert seq.translate() == "MA"

def test_translate_lowercase():
    seq = Sequence(
        id="seq1",
        sequence="atggcc"
    )

    assert seq.translate() == "MA"

def test_translate_iupac_unambiguous():
    seq = Sequence(
        id="seq1",
        sequence="ATGGCN"
    )

    assert seq.translate() == "MA"

def test_expand_iupac():
    assert expand_iupac("GCN") == ["GCA", "GCC", "GCG", "GCT"]

def test_translate_iupac_ambiguous():
    seq = Sequence(
        id="seq1",
        sequence="ATGRNY"
    )

    assert seq.translate() == "MX"