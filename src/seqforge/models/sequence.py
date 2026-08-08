from dataclasses import dataclass

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
            if letter in 'Aa':
                rev_compl += 'T'
            elif letter in 'Cc':
                rev_compl += 'G'
            elif letter in 'Gg':
                rev_compl += 'C'
            elif letter in 'Tt':
                rev_compl += 'A'

        return rev_compl
