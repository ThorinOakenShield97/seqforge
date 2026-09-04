from seqforge.models.sequence import Sequence
from seqforge.models.molecule_type import MoleculeType
from seqforge.models.filters import filter_by_length, filter_by_motif,filter_by_quality
from seqforge.models.fastq_read import FastqRead

import pytest

def test_filter_by_min_length():
    sequences = [
        Sequence(id="short", sequence="AT"),
        Sequence(id="medium", sequence="ATGC"),
        Sequence(id="long", sequence="ATGCGG"),
    ]

    result = filter_by_length(sequences, min_length=4)

    assert [sequence.id for sequence in result] == ["medium", "long"]


def test_filter_by_max_length():
    sequences = [
        Sequence(id="short", sequence="AT"),
        Sequence(id="medium", sequence="ATGC"),
        Sequence(id="long", sequence="ATGCGG"),
    ]

    result = filter_by_length(sequences, max_length=4)

    assert [sequence.id for sequence in result] == ["short", "medium"]


def test_filter_by_length_range_is_inclusive():
    sequences = [
        Sequence(id="below", sequence="AT"),
        Sequence(id="min", sequence="ATG"),
        Sequence(id="middle", sequence="ATGC"),
        Sequence(id="max", sequence="ATGCG"),
        Sequence(id="above", sequence="ATGCGG"),
    ]

    result = filter_by_length(
        sequences,
        min_length=3,
        max_length=5,
    )

    assert [sequence.id for sequence in result] == [
        "min",
        "middle",
        "max",
    ]


def test_filter_by_length_returns_empty_when_nothing_matches():
    sequences = [
        Sequence(id="one", sequence="AT"),
        Sequence(id="two", sequence="ATGC"),
    ]

    result = filter_by_length(sequences, min_length=10)

    assert result == []


def test_filter_by_length_rejects_invalid_range():
    sequences = [
        Sequence(id="seq1", sequence="ATGC"),
    ]

    with pytest.raises(ValueError):
        filter_by_length(
            sequences,
            min_length=10,
            max_length=5,
        )


def test_filter_by_length_rejects_negative_limits():
    sequences = [
        Sequence(id="seq1", sequence="ATGC"),
    ]

    with pytest.raises(ValueError):
        filter_by_length(sequences, min_length=-1)

    with pytest.raises(ValueError):
        filter_by_length(sequences, max_length=-1)

def test_filter_by_length_requires_at_least_one_limit():
    sequences = [
        Sequence(id="seq1", sequence="ATGC"),
    ]

    with pytest.raises(ValueError):
        filter_by_length(sequences)

def test_filter_by_motif():
    sequences = [
        Sequence(id="seq1", sequence="ATGCGT"),
        Sequence(id="seq2", sequence="GGGCCC"),
        Sequence(id="seq3", sequence="TTATGC"),
    ]

    result = filter_by_motif(sequences, "ATG")

    assert [sequence.id for sequence in result] == ["seq1", "seq3"]


def test_filter_by_motif_returns_empty_when_no_sequence_matches():
    sequences = [
        Sequence(id="seq1", sequence="ATGC"),
        Sequence(id="seq2", sequence="GGCC"),
    ]

    result = filter_by_motif(sequences, "TTT")

    assert result == []


def test_filter_by_motif_keeps_original_sequence_objects():
    sequences = [
        Sequence(id="seq1", sequence="ATGCGT"),
        Sequence(id="seq2", sequence="GGGCCC"),
    ]

    result = filter_by_motif(sequences, "ATG")

    assert result[0] is sequences[0]


def test_filter_by_motif_works_with_rna_and_protein():
    sequences = [
        Sequence(id="rna", sequence="AUGCGC", molecule_type=MoleculeType.RNA),
        Sequence(id="protein", sequence="MKWVTF", molecule_type=MoleculeType.PROTEIN),
    ]

    assert filter_by_motif(sequences, "AUG") == [sequences[0]]
    assert filter_by_motif(sequences, "WV") == [sequences[1]]


def test_filter_by_motif_does_not_modify_input_collection():
    sequences = [
        Sequence(id="seq1", sequence="ATGCGT"),
        Sequence(id="seq2", sequence="GGGCCC"),
    ]

    original = sequences.copy()

    filter_by_motif(sequences, "ATG")

    assert sequences == original


def test_filter_by_quality():
    reads = [
        FastqRead(id="low", sequence="ATGC", quality="!!!!"),  # 0
        FastqRead(id="medium", sequence="ATGC", quality="IIII"),  # 40
        FastqRead(id="high", sequence="ATGC", quality="JJJJ"),  # 41
    ]

    result = filter_by_quality(reads, min_quality=40)

    assert [read.id for read in result] == ["medium", "high"]


def test_filter_by_quality_includes_threshold():
    reads = [
        FastqRead(id="read1", sequence="ATGC", quality="IIII"),
    ]

    result = filter_by_quality(reads, min_quality=40)

    assert result == reads


def test_filter_by_quality_returns_empty_when_nothing_matches():
    reads = [
        FastqRead(id="low", sequence="ATGC", quality="!!!!"),
        FastqRead(id="medium", sequence="ATGC", quality="####"),
    ]

    result = filter_by_quality(reads, min_quality=40)

    assert result == []


def test_filter_by_quality_keeps_original_read_objects():
    reads = [
        FastqRead(id="read1", sequence="ATGC", quality="IIII"),
    ]

    result = filter_by_quality(reads, min_quality=30)

    assert result[0] is reads[0]


def test_filter_by_quality_rejects_negative_threshold():
    reads = [
        FastqRead(id="read1", sequence="ATGC", quality="IIII"),
    ]

    with pytest.raises(ValueError):
        filter_by_quality(reads, min_quality=-1)

def test_filter_by_quality_accepts_zero_threshold():
    reads = [
        FastqRead(id="read1", sequence="ATGC", quality="!!!!"),
    ]

    result = filter_by_quality(reads, min_quality=0)

    assert result == reads

def test_filter_by_quality_uses_mean_quality():
    reads = [
        FastqRead(
            id="read1",
            sequence="ATGC",
            quality="I?II",
        ),
        FastqRead(
            id="read2",
            sequence="ATGC",
            quality="IIII",
        ),
    ]

    result = filter_by_quality(reads, min_quality=39)

    assert [read.id for read in result] == ["read2"]

def test_filter_by_quality_with_empty_collection():
    result = filter_by_quality([], min_quality=30)

    assert result == []

def test_filter_by_quality_rejects_none_min_quality():
    reads = [
        FastqRead(id="read1", sequence="ATGC", quality="IIII"),
    ]

    with pytest.raises(ValueError):
        filter_by_quality(reads, min_quality=None)

def test_filter_by_motif_rejects_empty_motif():
    sequences = [
        Sequence(id="seq1", sequence="ATGC"),
    ]

    with pytest.raises(ValueError):
        filter_by_motif(sequences, motif="")