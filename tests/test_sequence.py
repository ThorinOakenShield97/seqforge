import pytest

from seqforge.models.sequence import (Sequence, expand_iupac, interpret_codon, CodonResultKind)

def test_sequence_is_public_api():
    from seqforge import Sequence

    assert Sequence is not None

def test_expand_iupac_is_public_api():
    from seqforge import expand_iupac

    assert expand_iupac("GCN") == [
        "GCA",
        "GCC",
        "GCG",
        "GCT",
    ]


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

def test_base_counts():
    seq = Sequence(
        id="seq1",
        sequence="ATGC"
    )

    assert seq.base_counts() == {
        "A": 1,
        "C": 1,
        "G": 1,
        "T": 1
    }

def test_base_counts_iupac():
    seq = Sequence(
        id="seq1",
        sequence="ATGCNRY"
    )

    assert seq.base_counts() == {
        "A": 1,
        "T": 1,
        "G": 1,
        "C": 1,
        "N": 1,
        "R": 1,
        "Y": 1
    }

def test_base_counts_lowercase():
    seq = Sequence(
        id="seq1",
        sequence="atgcnry"
    )

    assert seq.base_counts() == {
        "A": 1,
        "T": 1,
        "G": 1,
        "C": 1,
        "N": 1,
        "R": 1,
        "Y": 1
    }

def test_base_frequencies():
    seq = Sequence(
        id="seq1",
        sequence="ATGCC"
    )

    assert seq.base_frequencies() == {
        "A": 20.0,
        "T": 20.0,
        "G": 20.0,
        "C": 40.0
    }

def test_base_frequencies_empty_sequence():
    seq = Sequence(
        id="seq1",
        sequence=""
    )

    with pytest.raises(ValueError):
        seq.base_frequencies()

def test_base_frequencies_sum_to_100():
    seq = Sequence(
        id="seq1",
        sequence="ATGCC"
    )

    frequencies = seq.base_frequencies()

    assert sum(frequencies.values()) == 100.0

def test_find_motif():
    seq = Sequence(
        id="seq1",
        sequence="ATGCCATG"
    )

    assert seq.find_motif("ATG") == [0, 5]

def test_find_motif_overlapping():
    seq = Sequence(
        id="seq1",
        sequence="AAAA"
    )

    assert seq.find_motif("AAA") == [0, 1]

def test_find_motif_lowercase():
    seq = Sequence(
        id="seq1",
        sequence="atgccatg"
    )

    assert seq.find_motif("atg") == [0, 5]

def test_find_motif_empty():
    seq = Sequence(
        id="seq1",
        sequence="ATGC"
    )

    with pytest.raises(ValueError):
        seq.find_motif("")

def test_find_motif_not_found():
    seq = Sequence(
        id="seq1",
        sequence="ATGC"
    )

    assert seq.find_motif("AAA") == []

def test_find_orfs():
    seq = Sequence(
        id="seq1",
        sequence="ATGAAATAG"
    )

    assert seq.find_orfs() == ["ATGAAATAG"]

def test_find_orfs_offset():
    seq = Sequence(
        id="seq1",
        sequence="CATGAAATAG"
    )

    assert seq.find_orfs() == ["ATGAAATAG"]

def test_find_orfs_without_stop():
    seq = Sequence(
        id="seq1",
        sequence="ATGAAA"
    )

    assert seq.find_orfs() == []

def test_find_orfs_immediate_stop():
    seq = Sequence(
        id="seq1",
        sequence="ATGTAA"
    )

    assert seq.find_orfs() == ["ATGTAA"]

def test_find_multiple_orfs():
    seq = Sequence(
        id="seq1",
        sequence="ATGAAATAGCCCATGTTTTAA"
    )

    assert seq.find_orfs() == [
        "ATGAAATAG",
        "ATGTTTTAA"
    ]

def test_find_orfs_ignores_out_of_frame_stop():
    seq = Sequence(
        id="seq1",
        sequence="ATGAAATTA"
    )

    assert seq.find_orfs() == []

def test_find_orfs_frame_two():
    seq = Sequence(
        id="seq1",
        sequence="CCATGAAATAG"
    )

    assert seq.find_orfs() == ["ATGAAATAG"]

def test_find_orfs_with_iupac():
    seq = Sequence(
        id="seq1",
        sequence="ATGGCNTAA"
    )

    assert seq.find_orfs() == ["ATGGCNTAA"]

def test_find_orfs_ambiguous_possible_stop():
    seq = Sequence(
        id="seq1",
        sequence="ATGTANTAA"
    )

    assert seq.find_orfs() == ["ATGTANTAA"]

def test_find_orfs_iupac_stop():
    seq = Sequence(
        id="seq1",
        sequence="ATGTAR"
    )

    assert seq.find_orfs() == ["ATGTAR"]

def test_find_orfs_iupac_possible_stop():
    seq = Sequence(
        id="seq1",
        sequence="ATGTGRTAA"
    )

    assert seq.find_orfs() == ["ATGTGRTAA"]

def test_find_overlapping_orfs():
    seq = Sequence(
        id="seq1",
        sequence="ATGATGTAA"
    )

    assert seq.find_orfs() == [
        "ATGATGTAA",
        "ATGTAA"
    ]
def test_expand_iupac_unambiguous():
    assert expand_iupac("ATG") == ["ATG"]

def test_find_orfs_reverse_strand():
    seq = Sequence(
        id="seq1",
        sequence="TTACAT"
    )

    assert seq.find_orfs(strand="reverse") == ["ATGTAA"]


def test_find_orfs_invalid_strand():
    seq = Sequence(
        id="seq1",
        sequence="ATGTAA"
    )

    with pytest.raises(ValueError):
        seq.find_orfs(strand="banana")

def test_find_orfs_both_strands():
    seq = Sequence(
        id="seq1",
        sequence="TTACAT"
    )

    assert seq.find_orfs(strand="both") == ["ATGTAA"]

def test_find_orfs_both_strands_returns_forward_then_reverse():
    seq = Sequence(
        id="seq1",
        sequence="TTACATGTAG"
    )

    assert seq.find_orfs(strand="both") == [
        "ATGTAG",
        "ATGTAA"
    ]

def test_expand_iupac_rejects_empty_codon():
    with pytest.raises(ValueError):
        expand_iupac("")

def test_expand_iupac_rejects_incomplete_codon():
    with pytest.raises(ValueError):
        expand_iupac("AT")

def test_expand_iupac_rejects_long_codon():
    with pytest.raises(ValueError):
        expand_iupac("ATGG")

def test_expand_iupac_rejects_invalid_base():
    with pytest.raises(ValueError):
        expand_iupac("TXR")

def test_expand_iupac_accepts_lowercase():
    assert expand_iupac("tar") == ["TAA", "TAG"]

def test_translate_unambiguous_iupac_codon():
    seq = Sequence(
        id="seq1",
        sequence="ATGGCNTAA"
    )

    assert seq.translate() == "MA"

def test_translate_ambiguous_iupac_codon():
    seq = Sequence(
        id="seq1",
        sequence="ATGTANTAA"
    )

    assert seq.translate() == "MX"

def test_translate_ambiguous_stop_codon():
    seq = Sequence(
        id="seq1",
        sequence="ATGTARCCC"
    )

    assert seq.translate() == "M"

def test_translate_without_start_returns_empty_string():
    seq = Sequence(
        id="seq1",
        sequence="CCCGGGTTT"
    )

    assert seq.translate() == ""

def test_translate_frame_one():
    seq = Sequence(
        id="seq1",
        sequence="AATGAAATAG"
    )

    assert seq.translate(frame=1) == "MK"

def test_translate_without_frame_keeps_current_behavior():
    seq = Sequence(
        id="seq1",
        sequence="AATGAAATAG"
    )

    assert seq.translate() == "MK"

def test_translate_rejects_invalid_frame():
    seq = Sequence(
        id="seq1",
        sequence="ATGAAATAG"
    )

    with pytest.raises(ValueError):
        seq.translate(frame=3)

def test_translate_frame_two():
    seq = Sequence(
        id="seq1",
        sequence="CCATGAAATAG"
    )

    assert seq.translate(frame=2) == "MK"

def test_translate_frame_without_start_returns_empty_string():
    seq = Sequence(
        id="seq1",
        sequence="AATGAAATAG"
    )

    assert seq.translate(frame=0) == ""

def test_translate_none_frame_keeps_current_behavior():
    seq = Sequence(
        id="seq1",
        sequence="AATGAAATAG"
    )

    assert seq.translate(frame=None) == "MK"

def test_translate_uses_first_start_codon_in_frame():
    seq = Sequence(
        id="seq1",
        sequence="AATGAAATGCCCTAG"
    )

    assert seq.translate(frame=1) == "MKCP"

def test_find_orfs_frame_zero_without_start():
    seq = Sequence(
        id="seq1",
        sequence="AATGAAATAG"
    )

    assert seq.find_orfs(frame=0) == []

def test_find_orfs_frame_one():
    seq = Sequence(
        id="seq1",
        sequence="AATGAAATAG"
    )

    assert seq.find_orfs(frame=1) == ["ATGAAATAG"]

def test_find_orfs_both_strands_with_frame():
    seq = Sequence(
        id="seq1",
        sequence="AATGTAATTACATA"
    )

    assert seq.find_orfs(strand="both", frame=1) == [
        "ATGTAA",
        "ATGTAA"
    ]

def test_find_orfs_both_strands_with_frame_returns_distinct_orfs():
    seq = Sequence(
        id="seq1",
        sequence="CGACTACATGTGA"
    )

    assert seq.find_orfs(strand="both", frame=1) == [
        "ATGTGA",
        "ATGTAG",
    ]

def test_find_orfs_both_strands_respects_frame():
    seq = Sequence(
        id="seq1",
        sequence="CGACTACATGTGA"
    )

    assert seq.find_orfs(strand="both", frame=1) == [
        "ATGTGA",
        "ATGTAG",
    ]

def test_interpret_codon_returns_amino_acid_result():
    result = interpret_codon("GCN")

    assert result.kind == CodonResultKind.AMINO_ACID
    assert result.amino_acid == "A"

def test_interpret_codon_returns_stop_result():
    result = interpret_codon("TAR")

    assert result.kind == CodonResultKind.STOP
    assert result.amino_acid is None

def test_interpret_codon_returns_ambiguous_result():
    result = interpret_codon("TAN")

    assert result.kind == CodonResultKind.AMBIGUOUS
    assert result.amino_acid is None
