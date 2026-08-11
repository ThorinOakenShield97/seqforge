from typer.testing import CliRunner
from seqforge.commands.version import app

def test_version_command():
    runner = CliRunner()

    result = runner.invoke(app)

    assert result.exit_code == 0

from importlib.metadata import version

def test_version_command_shows_package_version():
    runner = CliRunner()

    result = runner.invoke(app)

    assert f"Version: {version('seqforge')}" in result.stdout