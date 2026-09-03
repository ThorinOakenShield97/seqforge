import pytest

from seqforge.models.sequence import (Sequence, expand_iupac, interpret_codon, CodonResultKind)
from seqforge.models.molecule_type import MoleculeType
from seqforge.models.filters import filter_by_length

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

    assert seq.translate(frame=2) == "MK"

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
        seq.translate(frame=0)

def test_translate_frame_two():
    seq = Sequence(
        id="seq1",
        sequence="CCATGAAATAG"
    )

    assert seq.translate(frame=3) == "MK"

def test_translate_frame_without_start_returns_empty_string():
    seq = Sequence(
        id="seq1",
        sequence="AATGAAATAG"
    )

    assert seq.translate(frame=1) == ""

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

    assert seq.translate(frame=2) == "MKCP"

def test_find_orfs_frame_zero_without_start():
    seq = Sequence(
        id="seq1",
        sequence="AATGAAATAG"
    )

    assert seq.find_orfs(frame=1) == []

def test_find_orfs_frame_one():
    seq = Sequence(
        id="seq1",
        sequence="AATGAAATAG"
    )

    assert seq.find_orfs(frame=2) == ["ATGAAATAG"]

def test_find_orfs_both_strands_with_frame():
    seq = Sequence(
        id="seq1",
        sequence="AATGTAATTACATA"
    )

    assert seq.find_orfs(strand="both", frame=2) == [
        "ATGTAA",
        "ATGTAA"
    ]

def test_find_orfs_both_strands_with_frame_returns_distinct_orfs():
    seq = Sequence(
        id="seq1",
        sequence="CGACTACATGTGA"
    )

    assert seq.find_orfs(strand="both", frame=2) == [
        "ATGTGA",
        "ATGTAG",
    ]

def test_find_orfs_both_strands_respects_frame():
    seq = Sequence(
        id="seq1",
        sequence="CGACTACATGTGA"
    )

    assert seq.find_orfs(strand="both", frame=2) == [
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

def test_transcribe_template_strand():
    seq = Sequence(
        id="seq1",
        sequence="GCAT",
    )

    assert seq.transcribe(strand="template") == "AUGC"

def test_transcribe_default_uses_coding_strand():
    seq = Sequence(
        id="seq1",
        sequence="ATGC",
    )

    assert seq.transcribe() == seq.transcribe(strand="coding")

def test_transcribe_rejects_invalid_strand():
    seq = Sequence(
        id="seq1",
        sequence="ATGC",
    )

    with pytest.raises(ValueError):
        seq.transcribe(strand="reverse")

def test_transcribe_template_strand_is_case_insensitive():
    seq = Sequence(
        id="seq1",
        sequence="gcat",
    )

    assert seq.transcribe(strand="template") == "AUGC"

def test_transcribe_template_matches_reverse_complement_transcription():
    seq = Sequence(
        id="seq1",
        sequence="GCAT",
    )

    assert seq.transcribe(strand="template") == (
        seq.reverse_complement().replace("T", "U")
    )

def test_translate_frame_three():
    seq = Sequence(
        id="seq1",
        sequence="CCATGAAATAG",
    )

    assert seq.translate(frame=3) == "MK"

def test_translate_rejects_frame_zero():
    seq = Sequence(
        id="seq1",
        sequence="ATGAAATAG",
    )

    with pytest.raises(ValueError):
        seq.translate(frame=0)

def test_find_orfs_frame_three():
    seq = Sequence(
        id="seq1",
        sequence="CCATGAAATAG",
    )

    assert seq.find_orfs(frame=3) == ["ATGAAATAG"]

def test_sequence_returns_kmers():
    seq = Sequence(
        id="seq1",
        sequence="ATGC",
    )

    assert seq.kmers(k=2) == [
        "AT",
        "TG",
        "GC",
    ]

def test_kmers_rejects_non_positive_k():
    seq = Sequence(
        id="seq1",
        sequence="ATGC",
    )

    with pytest.raises(ValueError):
        seq.kmers(k=0)

def test_kmers_rejects_negative_k():
    seq = Sequence(
        id="seq1",
        sequence="ATGC",
    )

    with pytest.raises(ValueError):
        seq.kmers(k=-1)

def test_kmers_returns_empty_for_k_larger_than_sequence():
    seq = Sequence(
        id="seq1",
        sequence="ATGC",
    )

    assert seq.kmers(k=5) == []

def test_kmer_counts():
    seq = Sequence(
        id="seq1",
        sequence="ATAT",
    )

    assert seq.kmer_counts(k=2) == {
        "AT": 2,
        "TA": 1,
    }

def test_kmer_counts_rejects_non_positive_k():
    seq = Sequence(
        id="seq1",
        sequence="ATGC",
    )

    with pytest.raises(ValueError):
        seq.kmer_counts(k=0)

def test_kmer_counts_returns_empty_for_k_larger_than_sequence():
    seq = Sequence(
        id="seq1",
        sequence="ATGC",
    )

    assert seq.kmer_counts(k=5) == {}

def test_sequence_can_be_rna():
    seq = Sequence(
        id="rna1",
        sequence="AUGC",
        molecule_type=MoleculeType.RNA,
    )

    assert seq.molecule_type == MoleculeType.RNA

def test_sequence_defaults_to_dna():
    seq = Sequence(
        id="seq1",
        sequence="ATGC",
    )

    assert seq.molecule_type == MoleculeType.DNA

def test_rna_sequence_accepts_rna_bases():
    seq = Sequence(
        id="rna1",
        sequence="AUGC",
        molecule_type=MoleculeType.RNA,
    )

    assert seq.sequence == "AUGC"

def test_rna_sequence_rejects_dna_thymine():
    with pytest.raises(ValueError):
        Sequence(
            id="rna1",
            sequence="ATGC",
            molecule_type=MoleculeType.RNA,
        )

def test_dna_sequence_rejects_rna_uracil():
    with pytest.raises(ValueError):
        Sequence(
            id="dna1",
            sequence="AUGC",
            molecule_type=MoleculeType.DNA,
        )

def test_rna_sequence_rejects_lowercase_thymine():
    with pytest.raises(ValueError):
        Sequence(
            id="rna1",
            sequence="augt",
            molecule_type=MoleculeType.RNA,
        )

def test_dna_sequence_rejects_lowercase_uracil():
    with pytest.raises(ValueError):
        Sequence(
            id="dna1",
            sequence="augc",
            molecule_type=MoleculeType.DNA,
        )

def test_protein_sequence_accepts_amino_acids():
    seq = Sequence(
        id="protein1",
        sequence="MKWVTF",
        molecule_type=MoleculeType.PROTEIN,
    )

    assert seq.sequence == "MKWVTF"

def test_protein_sequence_rejects_invalid_amino_acid():
    with pytest.raises(ValueError):
        Sequence(
            id="protein1",
            sequence="MKWVTFB",
            molecule_type=MoleculeType.PROTEIN,
        )

def test_protein_sequence_accepts_lowercase_amino_acids():
    seq = Sequence(
        id="protein1",
        sequence="mkwvtf",
        molecule_type=MoleculeType.PROTEIN,
    )

    assert seq.sequence == "mkwvtf"

def test_protein_sequence_rejects_reverse_complement():
    seq = Sequence(
        id="protein1",
        sequence="MKWVTF",
        molecule_type=MoleculeType.PROTEIN,
    )

    with pytest.raises(ValueError):
        seq.reverse_complement()

def test_protein_sequence_rejects_transcription():
    seq = Sequence(
        id="protein1",
        sequence="MKWVTF",
        molecule_type=MoleculeType.PROTEIN,
    )

    with pytest.raises(ValueError):
        seq.transcribe()

def test_rna_sequence_supports_reverse_complement():
    seq = Sequence(
        id="rna1",
        sequence="AUGC",
        molecule_type=MoleculeType.RNA,
    )

    assert seq.reverse_complement() == "GCAU"

def test_rna_reverse_complement_uses_uracil():
    seq = Sequence(
        id="rna1",
        sequence="AUGU",
        molecule_type=MoleculeType.RNA,
    )

    assert seq.reverse_complement() == "ACAU"

def test_rna_sequence_can_be_translated():
    seq = Sequence(
        id="rna1",
        sequence="AUGAAAUAG",
        molecule_type=MoleculeType.RNA,
    )

    assert seq.translate() == "MK"

def test_rna_sequence_stops_at_rna_stop_codon():
    seq = Sequence(
        id="rna1",
        sequence="AUGAAAUAA",
        molecule_type=MoleculeType.RNA,
    )

    assert seq.translate() == "MK"

def test_rna_sequence_can_be_translated_in_frame_two():
    seq = Sequence(
        id="rna1",
        sequence="AAUGAAAUAG",
        molecule_type=MoleculeType.RNA,
    )

    assert seq.translate(frame=2) == "MK"

def test_rna_sequence_can_be_translated_in_frame_three():
    seq = Sequence(
        id="rna1",
        sequence="AA AUG AAA UAG".replace(" ", ""),
        molecule_type=MoleculeType.RNA,
    )

    assert seq.translate(frame=3) == "MK"

def test_rna_sequence_rejects_transcription():
    seq = Sequence(
        id="rna1",
        sequence="AUGC",
        molecule_type=MoleculeType.RNA,
    )

    with pytest.raises(ValueError):
        seq.transcribe()

def test_protein_sequence_rejects_translation():
    seq = Sequence(
        id="protein1",
        sequence="MKWVTF",
        molecule_type=MoleculeType.PROTEIN,
    )

    with pytest.raises(ValueError):
        seq.translate()

def test_protein_sequence_rejects_gc_content():
    seq = Sequence(
        id="protein1",
        sequence="MKWVTF",
        molecule_type=MoleculeType.PROTEIN,
    )

    with pytest.raises(ValueError):
        seq.gc_content()

def test_protein_sequence_supports_length():
    seq = Sequence(
        id="protein1",
        sequence="MKWVTF",
        molecule_type=MoleculeType.PROTEIN,
    )

    assert seq.length() == 6

def test_protein_sequence_rejects_base_counts():
    seq = Sequence(
        id="protein1",
        sequence="MKWVTF",
        molecule_type=MoleculeType.PROTEIN,
    )

    with pytest.raises(ValueError):
        seq.base_counts()

def test_protein_sequence_rejects_base_frequencies():
    seq = Sequence(
        id="protein1",
        sequence="MKWVTF",
        molecule_type=MoleculeType.PROTEIN,
    )

    with pytest.raises(ValueError):
        seq.base_frequencies()

def test_protein_sequence_supports_motif_search():
    seq = Sequence(
        id="protein1",
        sequence="MKWVTF",
        molecule_type=MoleculeType.PROTEIN,
    )

    assert seq.find_motif("WV") == [2]

def test_rna_sequence_supports_gc_content():
    seq = Sequence(
        id="rna1",
        sequence="AUGC",
        molecule_type=MoleculeType.RNA,
    )

    assert seq.gc_content() == 50.0

def test_rna_sequence_supports_base_frequencies():
    seq = Sequence(
        id="rna1",
        sequence="AUGC",
        molecule_type=MoleculeType.RNA,
    )

    assert seq.base_frequencies() == {
        "A": 25.0,
        "C": 25.0,
        "G": 25.0,
        "U": 25.0,
    }

def test_rna_sequence_can_find_orfs():
    seq = Sequence(
        id="rna1",
        sequence="AUGAAAUAA",
        molecule_type=MoleculeType.RNA,
    )

    assert seq.find_orfs() == ["AUGAAAUAA"]

def test_rna_orf_preserves_uracil():
    seq = Sequence(
        id="rna1",
        sequence="AUGCCCUAG",
        molecule_type=MoleculeType.RNA,
    )

    assert seq.find_orfs() == ["AUGCCCUAG"]

def test_rna_sequence_can_find_orfs_on_reverse_strand():
    seq = Sequence(
        id="rna1",
        sequence="UUAUUUCAU",
        molecule_type=MoleculeType.RNA,
    )

    assert seq.find_orfs(strand="reverse") == ["AUGAAAUAA"]

def test_rna_sequence_can_find_orfs_on_both_strands():
    seq = Sequence(
        id="rna1",
        sequence="UUAUUUCAU",
        molecule_type=MoleculeType.RNA,
    )

    assert seq.find_orfs(strand="both") == [
        "AUGAAAUAA",
    ]

def test_protein_sequence_counts_amino_acids():
    seq = Sequence(
        id="protein1",
        sequence="MKWVTFM",
        molecule_type=MoleculeType.PROTEIN,
    )

    assert seq.amino_acid_counts() == {
        "M": 2,
        "K": 1,
        "W": 1,
        "V": 1,
        "T": 1,
        "F": 1,
    }

def test_protein_sequence_calculates_amino_acid_frequencies():
    seq = Sequence(
        id="protein1",
        sequence="MKWVTFM",
        molecule_type=MoleculeType.PROTEIN,
    )

    frequencies = seq.amino_acid_frequencies()

    assert frequencies["M"] == 2 / 7 * 100
    assert frequencies["K"] == 1 / 7 * 100
    assert frequencies["W"] == 1 / 7 * 100
    assert frequencies["V"] == 1 / 7 * 100
    assert frequencies["T"] == 1 / 7 * 100
    assert frequencies["F"] == 1 / 7 * 100

def test_protein_sequence_amino_acid_frequency_for_repeated_residue():
    seq = Sequence(
        id="protein1",
        sequence="MMMM",
        molecule_type=MoleculeType.PROTEIN,
    )

    assert seq.amino_acid_frequencies() == {
        "M": 100.0,
    }

def test_protein_sequence_supports_kmers():
    seq = Sequence(
        id="protein1",
        sequence="MKWVTF",
        molecule_type=MoleculeType.PROTEIN,
    )

    assert seq.kmers(k=2) == [
        "MK",
        "KW",
        "WV",
        "VT",
        "TF",
    ]

def test_dna_sequence_rejects_amino_acid_counts():
    seq = Sequence(
        id="dna1",
        sequence="ATGC",
        molecule_type=MoleculeType.DNA,
    )

    with pytest.raises(ValueError):
        seq.amino_acid_counts()

def test_protein_sequence_rejects_orf_detection():
    seq = Sequence(
        id="protein1",
        sequence="MKWVTF",
        molecule_type=MoleculeType.PROTEIN,
    )

    with pytest.raises(ValueError):
        seq.find_orfs()

def test_find_orfs_rna_supports_all_three_frames():
    seq = Sequence(
        id="rna1",
        sequence="AUGAAAUAGCCCAUGA",
        molecule_type=MoleculeType.RNA,
    )

    assert seq.find_orfs() == [
        "AUGAAAUAG",
    ]

def test_find_orfs_rna_supports_all_three_reading_frames():
    seq = Sequence(
        id="rna_frames",
        sequence="AAUAUGCCCUAAGCCAUGAAAUAG",
        molecule_type=MoleculeType.RNA,
    )

    orfs = seq.find_orfs()

    assert "AUGCCCUAA" in orfs
    assert "AUGAAAUAG" in orfs

def test_find_orfs_rna_frame_1():
    seq = Sequence(
        id="rna_f1",
        sequence="AUGCCCUAA",
        molecule_type=MoleculeType.RNA,
    )

    assert seq.find_orfs() == ["AUGCCCUAA"]


def test_find_orfs_rna_frame_2():
    seq = Sequence(
        id="rna_f2",
        sequence="AAUGCCCUAA",
        molecule_type=MoleculeType.RNA,
    )

    assert seq.find_orfs() == ["AUGCCCUAA"]


def test_find_orfs_rna_frame_3():
    seq = Sequence(
        id="rna_f3",
        sequence="AAAUGCCCUAA",
        molecule_type=MoleculeType.RNA,
    )

    assert seq.find_orfs() == ["AUGCCCUAA"]

def test_reverse_complement_rna():
    seq = Sequence(
        id="rna1",
        sequence="AUGC",
        molecule_type=MoleculeType.RNA,
    )

    assert seq.reverse_complement() == "GCAU"

def test_reverse_complement_protein_rejects_operation():
    seq = Sequence(
        id="protein1",
        sequence="MKWVTF",
        molecule_type=MoleculeType.PROTEIN,
    )

    with pytest.raises(ValueError):
        seq.reverse_complement()

def test_reverse_complement_dna():
    seq = Sequence(
        id="dna1",
        sequence="ATGC",
        molecule_type=MoleculeType.DNA,
    )

    assert seq.reverse_complement() == "GCAT"

def test_transcribe_rna_rejects_operation():
    seq = Sequence(
        id="rna1",
        sequence="AUGC",
        molecule_type=MoleculeType.RNA,
    )

    with pytest.raises(ValueError):
        seq.transcribe()

def test_transcribe_protein_rejects_operation():
    seq = Sequence(
        id="protein1",
        sequence="MKWVTF",
        molecule_type=MoleculeType.PROTEIN,
    )

    with pytest.raises(ValueError):
        seq.transcribe()

def test_transcribe_dna():
    seq = Sequence(
        id="dna1",
        sequence="ATGC",
        molecule_type=MoleculeType.DNA,
    )

    assert seq.transcribe() == "AUGC"

def test_translate_rna():
    seq = Sequence(
        id="rna1",
        sequence="AUGAAAUAG",
        molecule_type=MoleculeType.RNA,
    )

    assert seq.translate() == "MK"

def test_translate_protein_rejects_operation():
    seq = Sequence(
        id="protein1",
        sequence="MKWVTF",
        molecule_type=MoleculeType.PROTEIN,
    )

    with pytest.raises(ValueError):
        seq.translate()

def test_base_counts_rna():
    seq = Sequence(
        id="rna1",
        sequence="AUGCUU",
        molecule_type=MoleculeType.RNA,
    )

    assert seq.base_counts() == {
        "A": 1,
        "U": 3,
        "C": 1,
        "G": 1,
    }

def test_base_counts_protein_rejects_operation():
    seq = Sequence(
        id="protein1",
        sequence="MKWVTF",
        molecule_type=MoleculeType.PROTEIN,
    )

    with pytest.raises(ValueError):
        seq.base_counts()

def test_base_frequencies_rna():
    seq = Sequence(
        id="rna1",
        sequence="AUGCUU",
        molecule_type=MoleculeType.RNA,
    )

    assert seq.base_frequencies() == {
        "A": 1 / 6 * 100,
        "U": 3 / 6 * 100,
        "C": 1 / 6 * 100,
        "G": 1 / 6 * 100,
    }

def test_base_frequencies_protein_rejects_operation():
    seq = Sequence(
        id="protein1",
        sequence="MKWVTF",
        molecule_type=MoleculeType.PROTEIN,
    )

    with pytest.raises(ValueError):
        seq.base_frequencies()

def test_gc_content_rna():
    seq = Sequence(
        id="rna1",
        sequence="AUGCGC",
        molecule_type=MoleculeType.RNA,
    )

    assert seq.gc_content() == pytest.approx(66.67, rel=1e-4)

def test_gc_content_protein_rejects_operation():
    seq = Sequence(
        id="protein1",
        sequence="MKWVTF",
        molecule_type=MoleculeType.PROTEIN,
    )

    with pytest.raises(ValueError):
        seq.gc_content()

def test_length_protein():
    seq = Sequence(
        id="protein1",
        sequence="MKWVTF",
        molecule_type=MoleculeType.PROTEIN,
    )

    assert seq.length() == 6

def test_rna_sequence_rejects_thymine():
    with pytest.raises(ValueError):
        Sequence(
            id="rna_invalid",
            sequence="AUGCTU",
            molecule_type=MoleculeType.RNA,
        )

def test_dna_sequence_rejects_uracil():
    with pytest.raises(ValueError):
        Sequence(
            id="dna_invalid",
            sequence="ATGCU",
            molecule_type=MoleculeType.DNA,
        )

def test_rna_sequence_rejects_amino_acid_counts():
    seq = Sequence(
        id="rna1",
        sequence="AUGC",
        molecule_type=MoleculeType.RNA,
    )

    with pytest.raises(ValueError):
        seq.amino_acid_counts()

def test_rna_sequence_rejects_amino_acid_frequencies():
    seq = Sequence(
        id="rna1",
        sequence="AUGC",
        molecule_type=MoleculeType.RNA,
    )

    with pytest.raises(ValueError):
        seq.amino_acid_frequencies()

def test_sequence_accepts_string_molecule_type():
    seq = Sequence(
        id="dna1",
        sequence="ATGC",
        molecule_type="dna",
    )

    assert seq.molecule_type == MoleculeType.DNA

def test_sequence_rejects_invalid_molecule_type():
    with pytest.raises(ValueError):
        Sequence(
            id="invalid",
            sequence="ATGC",
            molecule_type="pepito",
        )

def test_sequence_accepts_uppercase_string_molecule_type():
    seq = Sequence(
        id="dna1",
        sequence="ATGC",
        molecule_type="DNA",
    )

    assert seq.molecule_type == MoleculeType.DNA

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