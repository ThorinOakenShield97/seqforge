from dataclasses import dataclass
from enum import Enum
from pathlib import Path

from seqforge.models.sequence import Sequence, MoleculeType
from seqforge.models.fastq_read import FastqRead
from seqforge.parsers.fasta import parse_fasta
from seqforge.parsers.fastq import parse_fastq

FASTA_EXTENSIONS = {".fasta", ".fa", ".fna"}
FASTQ_EXTENSIONS = {".fastq", ".fq"}


class InputSource(Enum):
    LITERAL = "literal"
    FASTA = "fasta"
    FASTQ = "fastq"

@dataclass
class ResolvedInput:
    sequences: list[Sequence]
    source: InputSource
    records: list[Sequence | FastqRead]

def resolve_input(value: str, molecule_type: MoleculeType | str = MoleculeType.DNA) -> ResolvedInput:
    """Resolve a literal sequence or FASTA file into Sequence objects."""
    input_path = Path(value)

    if isinstance(molecule_type, str):
        molecule_type = molecule_type.lower()
        molecule_type = MoleculeType(molecule_type)
    elif isinstance(molecule_type, MoleculeType):
        pass
    else:
        raise ValueError('Incorrect molecule type')

    if input_path.is_file():
        if input_path.suffix.lower() in FASTQ_EXTENSIONS:
            records = parse_fastq(input_path)
            return ResolvedInput(sequences = [], source = InputSource.FASTQ, records = records)
        # FASTA support will be added next.
        records = parse_fasta(input_path, molecule_type)
        return ResolvedInput(sequences = parse_fasta(input_path, molecule_type), source = InputSource.FASTA, records = records)

    if input_path.suffix.lower() in FASTA_EXTENSIONS:
        raise FileNotFoundError(f"Input file not found: {value}")

    record = Sequence(id='cli', sequence = value, molecule_type = molecule_type)
    return ResolvedInput(sequences = [record], source = InputSource.LITERAL, records = [record])
