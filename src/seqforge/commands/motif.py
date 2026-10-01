from seqforge.commands.input import InputSource, resolve_input
from seqforge.models.sequence import Sequence
from seqforge.models.molecule_type import MoleculeType

def motif(sequence:str, pattern:str):

    results = resolve_input(sequence)

    if results.source == InputSource.FASTQ:
        records = [Sequence(id = record.id, sequence = record.sequence, molecule_type = MoleculeType.DNA) for record in results.records]
    
    else:
        records = results.records

    for record in records:
        motifs = record.find_motif(pattern)
        print(f"positions: {motifs}")
