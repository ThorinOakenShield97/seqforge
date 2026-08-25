from dataclasses import dataclass
from itertools import product
from enum import Enum


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


CODON_TABLE = {
    'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
    'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
    'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
    'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
    'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
    'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
    'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
    'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
    'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
    'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
    'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
    'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
    'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
    'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
    'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
    'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
}

IUPAC_BASES = {
    'A': ['A'],
    'C': ['C'],
    'G': ['G'],
    'T': ['T'],
    'R': ['A', 'G'],
    'Y': ['C', 'T'],
    'S': ['G', 'C'],
    'W': ['A', 'T'],
    'K': ['G', 'T'],
    'M': ['A', 'C'],
    'B': ['C', 'G', 'T'],
    'D': ['A', 'G', 'T'],
    'H': ['A', 'C', 'T'],
    'V': ['A', 'C', 'G'],
    'N': ['A', 'C', 'G', 'T'],
}

def expand_iupac(codon:str) -> list[str]:
    """Return all possible DNA codon expansions for an IUPAC codon.

    The codon is case-insensitive and must contain exactly three valid
    IUPAC nucleotide symbols.

    Args:
        codon: Three-character DNA codon using IUPAC nucleotide symbols.

    Returns:
        A list containing all possible unambiguous DNA codons.

    Raises:
        ValueError: If the codon does not contain exactly three characters
            or contains an unsupported IUPAC symbol.
    """
    codon = codon.upper()
    if len(codon) == 3:
        for letter in codon:
            if letter not in IUPAC_BASES:
                raise ValueError('Letter not in IUPAC_BASES')
    else:
        raise ValueError('Incorrect number of letters')

    possibilities = [IUPAC_BASES[letter] for letter in codon]
    return [''.join(combination) for combination in product(*possibilities)]


class CodonResultKind(Enum):
    AMINO_ACID = "amino_acid"
    STOP = "stop"
    AMBIGUOUS = "ambiguous"

@dataclass
class CodonResult:
    kind: CodonResultKind
    amino_acid: str | None = None

def interpret_codon(codon:str) -> CodonResult:
    results = expand_iupac(codon)
    aminoacids = [CODON_TABLE[codon] for codon in results]

    if all(aminoacid == '*' for aminoacid in aminoacids):
        return CodonResult(kind = CodonResultKind.STOP)

    if len(set(aminoacids)) == 1:
        return CodonResult(kind = CodonResultKind.AMINO_ACID, amino_acid = aminoacids[0])

    else:
        return CodonResult(kind = CodonResultKind.AMBIGUOUS)


@dataclass
class Sequence:
    id: str
    sequence: str

    def length(self) -> int:
        """Return the length of the sequence."""
        return len(self.sequence)

    def gc_content(self) -> float:
        """Return the GC content of the sequence as a percentage.

        Raises:
            ValueError: If the sequence is empty.
        """
        if not self.sequence:
            raise ValueError('Cannot calculate GC content of an empty sequence.')
        total = len(self.sequence)
        g_count = self.sequence.count('G') + self.sequence.count('g')
        c_count = self.sequence.count('C') + self.sequence.count('c')

        gc_count = g_count + c_count

        return gc_count/total * 100

    def base_counts(self) -> dict[str, int]:
        """Return the count of each symbol in the sequence, normalized to uppercase."""
        sequence = self.sequence.upper()
        counts = {}
        for letter in sequence:
            if letter in counts:
                counts[letter] += 1
            else:
                counts[letter] = 1
        return counts

    def base_frequencies(self) -> dict[str, float]:
        """Return the frequency of each symbol in the sequence as a percentage.

            Raises:
                ValueError: If the sequence is empty.
        """
        counts = self.base_counts()
        total = len(self.sequence)

        if total == 0:
            raise ValueError('Cannot calculate base frequencies of an empty sequence')

        frequencies = {}
        for base in counts:
            frequencies[base] = counts[base]/total * 100
        return frequencies

    def reverse_complement(self) -> str:
        """Return the reverse complement of the DNA sequence.

            Raises:
                ValueError: If the sequence contains an unsupported DNA base.
        """
        reverse = self.sequence[::-1]
        rev_compl = ''
        for letter in reverse:
            if letter.upper() in DNA_COMPLEMENTS:
                rev_compl += DNA_COMPLEMENTS[letter.upper()]
            else:
                raise ValueError("Invalid DNA base.")

        return rev_compl

    def transcribe(self, strand: str = 'coding' ) -> str:
        """Return the RNA transcript of the DNA sequence.

        Args:
            strand: DNA strand to transcribe: ``"coding"`` or ``"template"``.
                The coding strand is transcribed directly. The template strand
                is reverse-complemented before transcription.

        Returns:
            The RNA transcript.

        Raises:
            ValueError: If ``strand`` is not ``"coding"`` or ``"template"``.
        """
        if strand == 'coding':
            sequence = self.sequence.upper()
        elif strand == 'template':
            sequence = self.reverse_complement()
        else:
            raise ValueError('Invalid Strand')

        return sequence.replace('T', 'U')


    def translate(self, frame: int | None = None) -> str:     
        """
        Translate the sequence from the first start codon to the first stop codon.

        When ``frame`` is ``None``, the first start codon in the sequence is used.
        When ``frame`` is 0, 1, or 2, the first start codon in the selected
        reading frame is used.

        IUPAC-ambiguous codons are translated when all possible expansions
        produce the same amino acid. Ambiguous codons producing different
        results are represented by ``X``.

        Args:
            frame: Reading frame to search for a start codon, or ``None`` to
            search the entire sequence.

        Returns:
            The translated protein sequence, or an empty string if no start
            codon is found.

        Raises:
            ValueError: If ``frame`` is not ``None``, 0, 1, or 2.
         """
        
        sequence = self.sequence.upper()
        if frame is None:
            start = sequence.find('ATG')

        else:
            if frame > 2 or frame < 0:
                raise ValueError('Invalid Frame')
            start = -1
            for i in range(frame,len(sequence)-2,3):
                if sequence[i:i+3] =='ATG':
                    start = i
                    break

        if start != -1:
            protein = ''
            for i in range(start,len(sequence)-2,3):
                triplet = sequence[i:i+3]
                result = interpret_codon(triplet)
                if result.kind == CodonResultKind.STOP:
                    break
                elif result.kind == CodonResultKind.AMINO_ACID:
                    protein += result.amino_acid
                else:
                    protein += 'X'
        else:
            return ''
                
        return protein

    def find_motif(self,motif:str) -> list[int]:
        """Return the starting positions of all occurrences of a motif.
        
        The search is case-insensitive and includes overlapping matches.
        
        Raises:
            ValueError: If the motif is empty.
        """
        positions = []
        sequence = self.sequence.upper()
        motif = motif.upper()
        if len(motif) == 0:
            raise ValueError('Empty Motif')
        for i in range(len(sequence)):
            if sequence[i:i+len(motif)] == motif:
                positions.append(i)

        return positions

    def find_orfs(self, strand: str = "forward", frame: int | None = None) -> list[str]:
        """Return open reading frames found in the selected DNA strand.

        Searches for complete ORFs beginning with ``ATG`` and ending at a
        compatible stop codon. When ``frame`` is ``None``, all reading frames
        are searched. When ``frame`` is 0, 1, or 2, the search is restricted
        to that reading frame.

        IUPAC-ambiguous codons are supported and are considered stop codons
        only when all possible expansions are stop codons.

        Args:
            strand: Strand to search: ``"forward"``, ``"reverse"``, or
                ``"both"``. ``"both"`` returns forward ORFs followed by
                reverse ORFs.
            frame: Reading frame to search, or ``None`` to search all frames.
                The selected frame is applied independently to each strand.

        Returns:
            A list of complete ORF nucleotide sequences, including their
            start and stop codons.

        Raises:
            ValueError: If ``strand`` is not ``"forward"``, ``"reverse"``,
                or ``"both"``, or if ``frame`` is not ``None``, 0, 1, or 2.
        """

        results = []
        if strand == 'both':
            return self.find_orfs('forward',frame=frame) + self.find_orfs('reverse',frame=frame)
        if strand == "forward":
            sequence = self.sequence.upper()
        elif strand == "reverse":
            sequence = self.reverse_complement()
        else:
            raise ValueError('Invalid Strand')


        if frame is None:
            positions = range(len(sequence))
        elif 0 <= frame <= 2:
            positions = range(frame,len(sequence),3)
        else:
            raise ValueError('Invalid Frame')

        for i in positions:
            triplet = sequence[i:i+3]
            if triplet == 'ATG':
                orf = ''
                orf += triplet
                for j in range(i+3,len(sequence)-2,3):
                    triplet = sequence[j:j+3]
                    result = interpret_codon(triplet)
                    if result.kind == CodonResultKind.AMINO_ACID:
                        orf += triplet
                    elif result.kind == CodonResultKind.STOP:
                        orf += triplet
                        results.append(orf)
                        break
                    else:
                        orf += triplet

        return results
