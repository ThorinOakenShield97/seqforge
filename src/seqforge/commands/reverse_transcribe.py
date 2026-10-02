from seqforge.commands.input import resolve_input
from seqforge.models.sequence import Sequence
from seqforge.models.molecule_type import MoleculeType


def reverse_transcribe(sequence:str):
    
    results = resolve_input(sequence, molecule_type = MoleculeType.RNA)
    
    records = [Sequence(id = record.id, sequence = record.sequence, molecule_type = MoleculeType.RNA) for record in results.records]
        

    for record in records:
        reverse = record.reverse_transcribe()
        print(reverse.sequence)