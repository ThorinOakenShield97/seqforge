from pathlib import Path
from seqforge.methods.sequence import Sequence

def parse_fasta(path: Path) -> list[Sequence]:
    """Parse a FASTA file and return its sequences."""
    records = []
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
                current_seq += line.strip()

        # Save the last sequence after reaching the end of the file.
        my_seq = Sequence(id = current_id, sequence = current_seq)
        records.append(my_seq)

        
    return records

"""
1. Crear una lista vacía para guardar las secuencias.

2. Abrir el archivo.

3. Leer línea por línea.

4. Si la línea empieza por '>':
      - Si ya estábamos construyendo una secuencia...
            guardarla en la lista.
      - Empezar una nueva secuencia.

5. Si la línea no empieza por '>':
      - Añadir esa línea a la secuencia actual.

6. Cuando termine el archivo...
      - Guardar también la última secuencia.

7. Devolver la lista."""