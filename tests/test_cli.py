import pytest

from typer.testing import CliRunner
from seqforge.cli import app
from importlib.metadata import version
from seqforge.commands.input import resolve_input
from seqforge.exceptions import InvalidFastaError
from seqforge.commands.input import (InputSource, resolve_input)


def test_version_command():
    runner = CliRunner()

    result = runner.invoke(app, ["version"])

    assert result.exit_code == 0


def test_version_command_shows_package_version():
    runner = CliRunner()
    result = runner.invoke(app, ["version"])
    assert f"Version: {version('seqforge')}" in result.stdout


def test_gc_command():
    runner = CliRunner()
    result = runner.invoke(app, ["gc", "ATGC"])

    assert result.exit_code == 0
    assert result.stdout.strip() == "GC content: 50.0%"

def test_gc_command_rejects_empty_sequence():
    runner = CliRunner()

    result = runner.invoke(app, ["gc", ""])

    assert result.exit_code != 0
    assert "Cannot calculate GC content of an empty sequence." in result.output

def test_translate_command():
    runner = CliRunner()

    result = runner.invoke(app, ["translate", "ATGAAATAG"])

    assert result.exit_code == 0
    assert result.stdout.strip() == "Protein: MK"

def test_translate_command_without_start_codon():
    runner = CliRunner()

    result = runner.invoke(app, ["translate", "CCCGGG"])

    assert result.exit_code == 0
    assert result.stdout.strip() == "No start codon found."

def test_translate_command_with_frame():
    runner = CliRunner()

    result = runner.invoke(
        app,
        ["translate", "AATGAAATAG", "--frame", "1"],
    )

    assert result.exit_code == 0
    assert result.stdout.strip() == "Protein: MK"

def test_orf_command():
    runner = CliRunner()

    result = runner.invoke(
        app,
        ["orf", "ATGAAATAG"],
    )

    assert result.exit_code == 0
    assert result.stdout.strip() == "ATGAAATAG"

def test_orf_command_with_frame():
    runner = CliRunner()

    result = runner.invoke(
        app,
        ["orf", "AATGAAATAG", "--frame", "1"],
    )

    assert result.exit_code == 0
    assert result.stdout.strip() == "ATGAAATAG"

def test_orf_command_with_reverse_strand():
    runner = CliRunner()

    result = runner.invoke(
        app,
        ["orf", "CGACTACATGTGA", "--strand", "reverse"],
    )

    assert result.exit_code == 0
    assert result.stdout.strip() == "ATGTAG"

def test_orf_command_with_both_strands():
    runner = CliRunner()

    result = runner.invoke(
        app,
        ["orf", "CGACTACATGTGA", "--strand", "both"],
    )

    assert result.exit_code == 0
    assert result.stdout.strip().splitlines() == [
        "Forward:",
        "ATGTGA",
        "Reverse:",
        "ATGTAG",
    ]

def test_orf_command_with_both_strands_and_frame():
    runner = CliRunner()

    result = runner.invoke(
        app,
        ["orf", "CGACTACATGTGA", "--strand", "both", "--frame", "1"],
    )

    assert result.exit_code == 0
    assert result.stdout.strip().splitlines() == [
        "Forward:",
        "ATGTGA",
        "Reverse:",
        "ATGTAG",
    ]

def test_orf_command_rejects_invalid_frame():
    runner = CliRunner()

    result = runner.invoke(
        app,
        ["orf", "ATGAAATAG", "--frame", "3"],
    )

    assert result.exit_code != 0
    assert "Invalid Frame" in result.output

def test_stats_command():
    runner = CliRunner()

    result = runner.invoke(
        app,
        ["stats", "ATGC"],
    )

    assert result.exit_code == 0
    assert result.stdout.strip() == (
        "Length: 4\n"
        "GC content: 50.0%\n"
        "Base counts:\n"
        "A: 1\n"
        "C: 1\n"
        "G: 1\n"
        "T: 1\n"
        "Base frequencies:\n"
        "A: 25.0%\n"
        "C: 25.0%\n"
        "G: 25.0%\n"
        "T: 25.0%"
    )

def test_stats_command_rejects_empty_sequence():
    runner = CliRunner()

    result = runner.invoke(app, ["stats", ""])

    assert result.exit_code != 0
    assert "Cannot calculate statistics for an empty sequence." in result.output

def test_stats_command_accepts_fasta_file(tmp_path):
    fasta = tmp_path / "sequence.fasta"
    fasta.write_text(">seq1\nATGC\n")

    runner = CliRunner()
    result = runner.invoke(app, ["stats", str(fasta)])

    assert result.exit_code == 0
    assert result.stdout.strip() == (
        ">seq1\n"
        "Length: 4\n"
        "GC content: 50.0%\n"
        "Base counts:\n"
        "A: 1\n"
        "C: 1\n"
        "G: 1\n"
        "T: 1\n"
        "Base frequencies:\n"
        "A: 25.0%\n"
        "C: 25.0%\n"
        "G: 25.0%\n"
        "T: 25.0%"
    )

def test_stats_command_accepts_multiple_fasta_records(tmp_path):
    fasta = tmp_path / "sequences.fasta"
    fasta.write_text(
        ">seq1\n"
        "ATGC\n"
        ">seq2\n"
        "AATT\n"
    )

    runner = CliRunner()
    result = runner.invoke(app, ["stats", str(fasta)])

    assert result.exit_code == 0
    assert result.stdout.strip() == (
        ">seq1\n"
        "Length: 4\n"
        "GC content: 50.0%\n"
        "Base counts:\n"
        "A: 1\n"
        "C: 1\n"
        "G: 1\n"
        "T: 1\n"
        "Base frequencies:\n"
        "A: 25.0%\n"
        "C: 25.0%\n"
        "G: 25.0%\n"
        "T: 25.0%\n"
        ">seq2\n"
        "Length: 4\n"
        "GC content: 0.0%\n"
        "Base counts:\n"
        "A: 2\n"
        "T: 2\n"
        "Base frequencies:\n"
        "A: 50.0%\n"
        "T: 50.0%"
    )

def test_stats_command_rejects_missing_input_file():
    runner = CliRunner()

    result = runner.invoke(
        app,
        ["stats", "does-not-exist.fasta"],
    )

    assert result.exit_code != 0
    assert "Input file not found: does-not-exist.fasta" in result.output

def test_stats_command_rejects_invalid_fasta(tmp_path):
    fasta = tmp_path / "invalid.fasta"
    fasta.write_text("ATGC\n")

    runner = CliRunner()
    result = runner.invoke(app, ["stats", str(fasta)])

    assert result.exit_code != 0
    assert "Missing FASTA header." in result.output

def test_resolve_input_accepts_literal_sequence():
    resolved = resolve_input("ATGC")

    assert resolved.source == InputSource.LITERAL
    assert len(resolved.sequences) == 1
    assert resolved.sequences[0].id == "cli"
    assert resolved.sequences[0].sequence == "ATGC"

def test_resolve_input_reads_fasta(tmp_path):
    fasta = tmp_path / "sequence.fasta"
    fasta.write_text(">seq1\nATGC\n")

    resolved = resolve_input(str(fasta))

    assert resolved.source == InputSource.FASTA
    assert len(resolved.sequences) == 1
    assert resolved.sequences[0].id == "seq1"
    assert resolved.sequences[0].sequence == "ATGC"

def test_resolve_input_rejects_missing_fasta_file():
    with pytest.raises(FileNotFoundError, match="Input file not found: missing.fasta"):
        resolve_input("missing.fasta")

def test_resolve_input_rejects_invalid_fasta(tmp_path):
    fasta = tmp_path / "invalid.fasta"
    fasta.write_text("ATGC\n")

    with pytest.raises(InvalidFastaError):
        resolve_input(str(fasta))

def test_gc_command_accepts_fasta_file(tmp_path):
    fasta = tmp_path / "sequence.fasta"
    fasta.write_text(">seq1\nATGC\n")

    runner = CliRunner()
    result = runner.invoke(app, ["gc", str(fasta)])

    assert result.exit_code == 0
    assert result.stdout.strip() == (
        ">seq1\n"
        "GC content: 50.0%"
    )

def test_gc_command_accepts_multiple_fasta_records(tmp_path):
    fasta = tmp_path / "sequences.fasta"
    fasta.write_text(
        ">seq1\n"
        "ATGC\n"
        ">seq2\n"
        "AATT\n"
    )

    runner = CliRunner()
    result = runner.invoke(app, ["gc", str(fasta)])

    assert result.exit_code == 0
    assert result.stdout.strip().splitlines() == [
        ">seq1",
        "GC content: 50.0%",
        ">seq2",
        "GC content: 0.0%",
    ]
def test_translate_command_accepts_fasta_file(tmp_path):
    fasta = tmp_path / "sequence.fasta"
    fasta.write_text(">gene1\nATGAAATAG\n")

    runner = CliRunner()
    result = runner.invoke(app, ["translate", str(fasta)])

    assert result.exit_code == 0
    assert result.stdout.strip() == (
        ">gene1\n"
        "Protein: MK"
    )

def test_translate_command_accepts_multiple_fasta_records(tmp_path):
    fasta = tmp_path / "sequences.fasta"
    fasta.write_text(
        ">gene1\n"
        "ATGAAATAG\n"
        ">gene2\n"
        "ATGCCCTAA\n"
    )

    runner = CliRunner()
    result = runner.invoke(app, ["translate", str(fasta)])

    assert result.exit_code == 0
    assert result.stdout.strip().splitlines() == [
        ">gene1",
        "Protein: MK",
        ">gene2",
        "Protein: MP",
    ]

def test_translate_command_handles_missing_start_in_fasta(tmp_path):
    fasta = tmp_path / "sequences.fasta"
    fasta.write_text(
        ">gene1\n"
        "ATGAAATAG\n"
        ">gene2\n"
        "CCCGGG\n"
    )

    runner = CliRunner()
    result = runner.invoke(app, ["translate", str(fasta)])

    assert result.exit_code == 0
    assert result.stdout.strip().splitlines() == [
        ">gene1",
        "Protein: MK",
        ">gene2",
        "No start codon found.",
    ]

def test_orf_command_accepts_fasta_file(tmp_path):
    fasta = tmp_path / "sequence.fasta"
    fasta.write_text(">seq1\nATGAAATAG\n")

    runner = CliRunner()
    result = runner.invoke(app, ["orf", str(fasta)])

    assert result.exit_code == 0
    assert result.stdout.strip() == (
        ">seq1\n"
        "ATGAAATAG"
    )

def test_orf_command_accepts_multiple_fasta_records(tmp_path):
    fasta = tmp_path / "sequences.fasta"
    fasta.write_text(
        ">seq1\n"
        "ATGAAATAG\n"
        ">seq2\n"
        "ATGCCCTAA\n"
    )

    runner = CliRunner()
    result = runner.invoke(app, ["orf", str(fasta)])

    assert result.exit_code == 0
    assert result.stdout.strip().splitlines() == [
        ">seq1",
        "ATGAAATAG",
        ">seq2",
        "ATGCCCTAA",
    ]

def test_orf_command_with_both_strands_labels_results():
    runner = CliRunner()

    result = runner.invoke(
        app,
        ["orf", "CGACTACATGTGA", "--strand", "both", "--frame", "1"],
    )

    assert result.exit_code == 0
    assert result.stdout.strip().splitlines() == [
        "Forward:",
        "ATGTGA",
        "Reverse:",
        "ATGTAG",
    ]

def test_orf_command_with_both_strands_without_reverse_orf():
    runner = CliRunner()

    result = runner.invoke(
        app,
        ["orf", "AATGTAAAA", "--strand", "both", "--frame", "1"],
    )

    assert result.exit_code == 0
    assert result.stdout.strip().splitlines() == [
        "Forward:",
        "ATGTAA",
        "Reverse:",
        "No ORFs found.",
    ]

def test_gc_command_rejects_invalid_fasta(tmp_path):
    fasta = tmp_path / "invalid.fasta"
    fasta.write_text("ATGC\n")

    runner = CliRunner()
    result = runner.invoke(app, ["gc", str(fasta)])

    assert result.exit_code != 0
    assert "Missing FASTA header." in result.output

def test_translate_command_rejects_invalid_fasta(tmp_path):
    fasta = tmp_path / "invalid.fasta"
    fasta.write_text("ATGC\n")

    runner = CliRunner()
    result = runner.invoke(app, ["translate", str(fasta)])

    assert result.exit_code != 0
    assert "Missing FASTA header." in result.output

def test_translate_command_rejects_missing_input_file():
    runner = CliRunner()

    result = runner.invoke(
        app,
        ["translate", "missing.fasta"],
    )

    assert result.exit_code != 0
    assert "Input file not found: missing.fasta" in result.output

def test_orf_command_rejects_invalid_fasta(tmp_path):
    fasta = tmp_path / "invalid.fasta"
    fasta.write_text("ATGC\n")

    runner = CliRunner()
    result = runner.invoke(app, ["orf", str(fasta)])

    assert result.exit_code != 0
    assert "Missing FASTA header." in result.output

def test_orf_command_rejects_missing_input_file():
    runner = CliRunner()

    result = runner.invoke(
        app,
        ["orf", "missing.fasta"],
    )

    assert result.exit_code != 0
    assert "Input file not found: missing.fasta" in result.output