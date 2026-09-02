from pathlib import Path
from seqforge.models.sequence import Sequence,MoleculeType
from seqforge.exceptions import InvalidFastaError

def parse_fasta(path: Path, molecule_type: MoleculeType = MoleculeType.DNA) -> list[Sequence]:
    """Parse a FASTA file and return its sequences.
    Args:
        path: Path to the FASTA file.
        molecule_type: Molecule type assigned to each parsed sequence. Defaults to DNA.

    Returns:
        A list of parsed Sequence objects

    Raises:
        InvalidFastaError: If the FASTA structure is invalid.
        ValueError: If a parsed sequence is incompatible with molecule_type.
        """
    records = []
    current_id = None
    with open(path,'r') as file:
        current_seq = ''
        for line in file: 
            # A new FASTA record begins. Save the previous sequence before starting a new one.
            if line.startswith('>'):
                if len(current_seq) != 0:
                    my_seq = Sequence(id = current_id, sequence = current_seq, molecule_type = molecule_type)
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
        my_seq = Sequence(id = current_id, sequence = current_seq, molecule_type = molecule_type)
        records.append(my_seq)

        
    return records
