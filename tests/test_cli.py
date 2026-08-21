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