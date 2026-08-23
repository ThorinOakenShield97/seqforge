import platform

from importlib.metadata import version as package_version


def version() -> None:
    """Display SeqForge version information."""
    print("SeqForge")
    print(f"Version: {package_version('seqforge')}")
    print(f"Python: {platform.python_version()}")