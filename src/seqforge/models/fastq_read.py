from dataclasses import dataclass

@dataclass
class FastqRead:
    id: str
    sequence: str
    quality: str

    def __post_init__(self) -> None:
        if len(self.sequence) != len(self.quality):
            raise ValueError("Sequence and quality must have the same length.")
        if not self.sequence:
            raise ValueError("Sequence cannot be empty.")
        if any(ord(char) < 33 for char in self.quality):
            raise ValueError("Invalid FASTQ quality character.")

    def quality_scores(self) -> list[int]:
        return [ord(character) - 33 for character in self.quality]

    def mean_quality(self):
        results = self.quality_scores()
        total = len(results)
        return sum(results)/total

    def min_quality(self):
        results = self.quality_scores()
        return min(results)

    def max_quality(self):
        results = self.quality_scores()
        return max(results)