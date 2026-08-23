from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from seqforge.models.sequence import Sequence
from seqforge.parsers.fasta import parse_fasta

FASTA_EXTENSIONS = {".fasta", ".fa", ".fna"}


class InputSource(Enum):
    LITERAL = "literal"
    FASTA = "fasta"

@dataclass
class ResolvedInput:
    sequences: list[Sequence]
    source: InputSource

def resolve_input(value: str) -> ResolvedInput:
    """Resolve a literal sequence or FASTA file into Sequence objects."""
    input_path = Path(value)

    if input_path.is_file():
        # FASTA support will be added next.
        return ResolvedInput(sequences = parse_fasta(input_path), source=InputSource.FASTA)

    if input_path.suffix.lower() in FASTA_EXTENSIONS:
        raise FileNotFoundError(f"Input file not found: {value}")

    return ResolvedInput(sequences = [Sequence(id = "cli", sequence = value)], source = InputSource.LITERAL)