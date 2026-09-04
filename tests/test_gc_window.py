import pytest

from seqforge.models.sequence import Sequence
from seqforge.models.gc_window import gc_content_windows
from seqforge.models.molecule_type import MoleculeType


def test_gc_content_windows():
    sequence = Sequence(id="seq1", sequence="ATGCATGC")

    result = gc_content_windows(sequence, window_size=4)

    assert result == [
        (1, 4, 50.0),
        (2, 5, 50.0),
        (3, 6, 50.0),
        (4, 7, 50.0),
        (5, 8, 50.0),
    ]


def test_gc_content_windows_uses_sliding_windows():
    sequence = Sequence(id="seq1", sequence="GCGAAA")

    result = gc_content_windows(sequence, window_size=4)

    assert result == [
        (1, 4, 75.0),
        (2, 5, 50.0),
        (3, 6, 25.0),
    ]


def test_gc_content_windows_window_size_one():
    sequence = Sequence(id="seq1", sequence="ACGT")

    result = gc_content_windows(sequence, window_size=1)

    assert result == [
        (1, 1, 0.0),
        (2, 2, 100.0),
        (3, 3, 100.0),
        (4, 4, 0.0),
    ]

def test_gc_content_windows_window_size_equals_sequence_length():
    sequence = Sequence(id="seq1", sequence="ATGC")

    result = gc_content_windows(sequence, window_size=4)

    assert result == [
        (1, 4, 50.0),
    ]


def test_gc_content_windows_returns_empty_when_window_is_larger():
    sequence = Sequence(id="seq1", sequence="ATGC")

    result = gc_content_windows(sequence, window_size=5)

    assert result == []


def test_gc_content_windows_rejects_non_positive_window_size():
    sequence = Sequence(id="seq1", sequence="ATGC")

    with pytest.raises(ValueError):
        gc_content_windows(sequence, window_size=0)

    with pytest.raises(ValueError):
        gc_content_windows(sequence, window_size=-1)


def test_gc_content_windows_rejects_protein():
    sequence = Sequence(
        id="protein",
        sequence="MKWVTF",
        molecule_type=MoleculeType.PROTEIN,
    )

    with pytest.raises(ValueError):
        gc_content_windows(sequence, window_size=3)