from seqforge.commands.input import InputSource, resolve_input
from seqforge.models.sequence import Sequence
from seqforge.models.alignment import multiple_alignment

def align(sequence:str, molecule_type: str = "DNA"):

    results = resolve_input(sequence, molecule_type)

    if results.source == InputSource.FASTQ:
        records = [Sequence(id = record.id, sequence = record.sequence, molecule_type = molecule_type) for record in results.records]

    else:
        records = results.records
    
    alignments = multiple_alignment(records)

    for alignment in alignments:
        print(alignment)