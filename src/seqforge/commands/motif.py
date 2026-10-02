from seqforge.commands.input import InputSource, resolve_input
from seqforge.models.sequence import Sequence


def motif(sequence:str, pattern:str, molecule_type:str = 'DNA'):

    results = resolve_input(sequence, molecule_type)

    if results.source == InputSource.FASTQ:
        records = [Sequence(id = record.id, sequence = record.sequence, molecule_type = molecule_type) for record in results.records]
    
    else:
        records = results.records

    for record in records:
        motifs = record.find_motif(pattern)
        print(f"positions: {motifs}")
