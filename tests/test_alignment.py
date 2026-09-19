from seqforge.models.sequence import Sequence
from seqforge.models.alignment import global_alignment, multiple_alignment
from seqforge.models.molecule_type import MoleculeType

import pytest


def test_global_alignment():
    seq1 = Sequence(id="seq1", sequence="ATGC")
    seq2 = Sequence(id="seq2", sequence="ATTC")

    result = global_alignment(seq1, seq2)

    assert result == ("ATGC", "ATTC")

def test_global_alignment_inserts_gap():
    seq1 = Sequence(id="seq1", sequence="ATGC")
    seq2 = Sequence(id="seq2", sequence="ATC")

    result = global_alignment(seq1, seq2)

    assert result == ("ATGC", "AT-C")

def test_global_alignment_inserts_gap_in_seq1():
    seq1 = Sequence(id="seq1", sequence="ATC")
    seq2 = Sequence(id="seq2", sequence="ATGC")

    result = global_alignment(seq1, seq2)

    assert result == ("AT-C", "ATGC")

def test_global_alignment_custom_scoring():
    seq1 = Sequence(id="seq1", sequence="ATGC")
    seq2 = Sequence(id="seq2", sequence="ATC")

    result = global_alignment(
        seq1,
        seq2,
        match=2,
        mismatch=-1,
        gap=-2,
    )

    assert result == ("ATGC", "AT-C")

def test_global_alignment_rejects_different_molecule_types():
    seq1 = Sequence(
        id="dna1",
        sequence="ATGC",
        molecule_type=MoleculeType.DNA,
    )
    seq2 = Sequence(
        id="rna1",
        sequence="AUGC",
        molecule_type=MoleculeType.RNA,
    )

    with pytest.raises(ValueError):
        global_alignment(seq1, seq2)

def test_global_alignment_protein():
    seq1 = Sequence(
        id="protein1",
        sequence="MKWVTF",
        molecule_type=MoleculeType.PROTEIN,
    )
    seq2 = Sequence(
        id="protein2",
        sequence="MKWETF",
        molecule_type=MoleculeType.PROTEIN,
    )

    result = global_alignment(seq1, seq2)

    assert result == ("MKWVTF", "MKWETF")

def test_global_alignment_preserves_remaining_sequence():
    seq1 = Sequence(id="seq1", sequence="ACG")
    seq2 = Sequence(id="seq2", sequence="ATCGG")

    result = global_alignment(seq1, seq2)

    assert result == ("A-CG-", "ATCGG")


def test_multiple_alignment_two_sequences():
    seq1 = Sequence(id="seq1", sequence="ATGC")
    seq2 = Sequence(id="seq2", sequence="ATC")

    result = multiple_alignment([seq1, seq2])

    assert result == [
        "ATGC",
        "AT-C",
    ]

def test_multiple_alignment_three_sequences():
    seq1 = Sequence(id="seq1", sequence="ATGC")
    seq2 = Sequence(id="seq2", sequence="ATC")
    seq3 = Sequence(id="seq3", sequence="ATGC")

    result = multiple_alignment([seq1, seq2, seq3])

    assert result == [
        "ATGC",
        "AT-C",
        "ATGC",
    ]

def test_multiple_alignment_propagates_gaps():
    seq1 = Sequence(id="seq1", sequence="ATC")
    seq2 = Sequence(id="seq2", sequence="AC")
    seq3 = Sequence(id="seq3", sequence="AGTC")

    result = multiple_alignment([seq1, seq2, seq3])

    assert result == [
        "A-TC",
        "A--C",
        "AGTC",
    ]

def test_multiple_alignment_propagates_multiple_gaps():
    seq1 = Sequence(id="seq1", sequence="ACG", molecule_type=MoleculeType.DNA)
    seq2 = Sequence(id="seq2", sequence="ACG", molecule_type=MoleculeType.DNA)
    seq3 = Sequence(id="seq3", sequence="ATCGG", molecule_type=MoleculeType.DNA)

    result = multiple_alignment([seq1, seq2, seq3])

    assert result == [
        "A-CG-",
        "A-CG-",
        "ATCGG",
    ]

def test_multiple_alignment_four_sequences():
    seq1 = Sequence(id="seq1", sequence="ATC")
    seq2 = Sequence(id="seq2", sequence="AC")
    seq3 = Sequence(id="seq3", sequence="AGTC")
    seq4 = Sequence(id="seq4", sequence="ATGC")

    result = multiple_alignment([seq1, seq2, seq3, seq4])
    print(result)

    assert result == [
    "A-T-C",
    "A---C",
    "AGT-C",
    "A-TGC",
]

def test_multiple_alignment_single_sequence():
    seq1 = Sequence(id="seq1", sequence="ATGC")

    result = multiple_alignment([seq1])

    assert result == ["ATGC"]

def test_multiple_alignment_rejects_empty_list():
    with pytest.raises(ValueError):
        multiple_alignment([])

def test_multiple_alignment_rejects_different_molecule_types():
    seq1 = Sequence(
        id="seq1",
        sequence="ATGC",
        molecule_type=MoleculeType.DNA,
    )
    seq2 = Sequence(
        id="seq2",
        sequence="AUGC",
        molecule_type=MoleculeType.RNA,
    )

    with pytest.raises(ValueError):
        multiple_alignment([seq1, seq2])

def test_multiple_alignment_different_lengths():
    seq1 = Sequence(id="seq1", sequence="ACGT")
    seq2 = Sequence(id="seq2", sequence="AC")
    seq3 = Sequence(id="seq3", sequence="ACGTT")

    result = multiple_alignment([seq1, seq2, seq3])

    assert result == [
        "ACGT-",
        "AC---",
        "ACGTT",
    ]