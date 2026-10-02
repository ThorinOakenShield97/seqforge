from seqforge.models.sequence import Sequence

def distance(seq_1:str, seq_2:str, molecule_type:str = 'DNA'):

    seq_1 = Sequence(id = 'seq 1', sequence = seq_1, molecule_type = molecule_type)
    seq_2 = Sequence(id = 'seq 2', sequence = seq_2, molecule_type = molecule_type)

    result = seq_1.distance(seq_2)

    print(f"distance: {result}")