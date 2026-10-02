import pytest

from seqforge.models.sequence import Sequence
from seqforge.models.molecule_type import MoleculeType
from seqforge.models.distances import pairwise_distances


def test_pairwise_distances():
    seq1 = Sequence(id="seq1", sequence="AAAA")
    seq2 = Sequence(id="seq2", sequence="AAAT")
    seq3 = Sequence(id="seq3", sequence="AATT")

    result = pairwise_distances([seq1, seq2, seq3])

    assert result == [
        ("seq1", "seq2", 1),
        ("seq1", "seq3", 2),
        ("seq2", "seq3", 1),
    ]

def test_pairwise_distances_rejects_mixed_molecule_types():
    dna = Sequence(id="dna", sequence="ATGC", molecule_type=MoleculeType.DNA)
    rna = Sequence(id="rna", sequence="AUGC", molecule_type=MoleculeType.RNA)

    with pytest.raises(ValueError, match="same molecule type"):
        pairwise_distances([dna, rna])