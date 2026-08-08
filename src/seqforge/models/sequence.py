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
