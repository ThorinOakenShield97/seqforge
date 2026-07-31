from pathlib import Path
from seqforge.models.sequence import Sequence
from seqforge.exceptions import InvalidFastaError

def parse_fasta(path: Path) -> list[Sequence]:
    """Parse a FASTA file and return its sequences."""
    records = []
    current_id = None
    with open(path,'r') as file:
        current_seq = ''
        for line in file: 
            # A new FASTA record begins. Save the previous sequence before starting a new one.
            if line.startswith('>'):
                if len(current_seq) != 0:
                    my_seq = Sequence(id = current_id, sequence = current_seq)
                    records.append(my_seq)
                    current_seq = ''
                current_id = line[1:].strip()
            else:
                if current_id:
                    current_seq += line.strip()
                else:
                    raise InvalidFastaError('Missing FASTA header.')
 

        # Save the last sequence after reaching the end of the file.
        if len(current_seq) == 0:
            raise InvalidFastaError('Missing sequence')
        my_seq = Sequence(id = current_id, sequence = current_seq)
        records.append(my_seq)

        
    return records
