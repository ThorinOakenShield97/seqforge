from seqforge.commands.input import InputSource, resolve_input
from seqforge.models.molecule_type import MoleculeType
from seqforge.models.sequence import Sequence
from seqforge.models.alignment import multiple_alignment

def align(sequence:str):

    results = resolve_input(sequence)

    if results.source == InputSource.FASTQ:
        records = [Sequence(id = record.id, sequence = record.sequence, molecule_type = MoleculeType.DNA) for record in results.records]

    else:
        records = results.records
    
    alignments = multiple_alignment(records)

    for alignment in alignments:
        print(alignment)