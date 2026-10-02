import typer

from seqforge.commands.input import resolve_input, InputSource
from seqforge.models.sequence import Sequence
from seqforge.models.molecule_type import MoleculeType
from seqforge.models.distances import pairwise_distances

def distances(sequence:str, molecule_type:str = 'DNA'):

    results = resolve_input(sequence, molecule_type = molecule_type)

    if results.source == InputSource.FASTQ:
            records = [Sequence(id = record.id, sequence = record.sequence, molecule_type = molecule_type) for record in results.records]

    else:
        records = results.records
    
    distances = pairwise_distances(records)

    for seq_1, seq_2, distance in distances:
         print(f"{seq_1}, {seq_2}: {distance}")