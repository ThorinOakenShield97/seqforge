from pathlib import Path

from seqforge.parsers.fasta import parse_fasta

records = parse_fasta(Path("fasta.txt"))

for record in records:
    print(record)