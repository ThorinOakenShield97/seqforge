from typer.testing import CliRunner
from seqforge.commands.version import app
from seqforge.cli import app

def test_version_command():
    runner = CliRunner()

    result = runner.invoke(app, ["version"])

    assert result.exit_code == 0

from importlib.metadata import version



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
        "ATGTGA",
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
        "ATGTGA",
        "ATGTAG",
    ]

def test_orf_command_rejects_invalid_frame():
    runner = CliRunner()

    result = runner.invoke(
        app,
        ["orf", "ATGAAATAG", "--frame", "3"],
    )

    assert result.exit_code != 0

def test_orf_command_rejects_invalid_frame():
    runner = CliRunner()

    result = runner.invoke(
        app,
        ["orf", "ATGAAATAG", "--frame", "3"],
    )

    assert result.exit_code != 0
    assert "Invalid Frame" in result.output