from dataclasses import dataclass

@dataclass
class Sequence:
    id: str
    sequence: str

    def length(self):
        return len(self.sequence)
