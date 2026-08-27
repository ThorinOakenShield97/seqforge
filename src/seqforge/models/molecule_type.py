from enum import Enum
from dataclasses import dataclass

class MoleculeType(Enum):
    DNA = 'dna'
    RNA = 'rna'
    PROTEIN = 'protein'

