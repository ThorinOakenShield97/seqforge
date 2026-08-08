from dataclasses import dataclass


DNA_COMPLEMENTS = {
    'A': 'T',
    'T': 'A',
    'C': 'G',
    'G': 'C',
    'R': 'Y',
    'Y': 'R',
    'S': 'S',
    'W': 'W',
    'K': 'M',
    'M': 'K',
    'B': 'V',
    'V': 'B',
    'D': 'H',
    'H': 'D',
    'N': 'N'
}


@dataclass
class Sequence:
    id: str
    sequence: str

    def length(self):
        return len(self.sequence)

    def gc_content(self):
        if not self.sequence:
            raise ValueError('Cannot calculate GC content of an empty sequence.')
        total = len(self.sequence)
        g_count = self.sequence.count('G') + self.sequence.count('g')
        c_count = self.sequence.count('C') + self.sequence.count('c')

        gc_count = g_count + c_count

        return gc_count/total * 100

    def reverse_complement(self):
        reverse = self.sequence[::-1]
        rev_compl = ''
        for letter in reverse:
            if letter.upper() in DNA_COMPLEMENTS:
                rev_compl += DNA_COMPLEMENTS[letter.upper()]
            else:
                raise ValueError("Invalid DNA base.")

        return rev_compl

    def transcribe(self):
        return self.sequence.upper().replace('T', 'U')
