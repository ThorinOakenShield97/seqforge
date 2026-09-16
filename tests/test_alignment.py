from seqforge.models.sequence import Sequence
from seqforge.models.alignment import global_alignment
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