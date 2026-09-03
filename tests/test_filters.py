from seqforge.models.sequence import Sequence
from seqforge.models.molecule_type import MoleculeType
from seqforge.models.filters import filter_by_length, filter_by_motif

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