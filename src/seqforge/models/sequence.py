from dataclasses import dataclass
from itertools import product
from enum import Enum
from seqforge.models.molecule_type import MoleculeType


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

RNA_COMPLEMENTS = {
    'A': 'U',
    'U': 'A',
    'G': 'C',
    'C': 'G',
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

AMINOACIDS = {'A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W','Y'}

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
    molecule_type: MoleculeType = MoleculeType.DNA

    def __post_init__(self):

        if isinstance(self.molecule_type, MoleculeType):

            if self.molecule_type == MoleculeType.RNA:
                if 'T' in self.sequence.upper():
                    raise ValueError('RNA must not have Ts')
            elif self.molecule_type == MoleculeType.DNA:
                if 'U' in self.sequence.upper():
                        raise ValueError('DNA must not have Us')
            elif self.molecule_type == MoleculeType.PROTEIN:
                for letter in self.sequence.upper():
                    if letter not in AMINOACIDS:
                        raise ValueError('Not Valid Aminoacid in Protein')
                    
        elif isinstance(self.molecule_type, str):
            self.molecule_type = self.molecule_type.lower()
            self.molecule_type = MoleculeType(self.molecule_type)
        else:
            raise ValueError('Invalid molecule type')


    def length(self) -> int:
        """Return the length of the sequence."""
        return len(self.sequence)

    def gc_content(self) -> float:
        """Return the GC content of the sequence as a percentage.

        Raises:
            ValueError: If the sequence is empty or is a protein.
        """
        if not self.sequence:
            raise ValueError('Cannot calculate GC content of an empty sequence.')
        elif self.molecule_type == MoleculeType.PROTEIN:
            raise ValueError('Proteins do not have Gs or Cs')
        total = len(self.sequence)
        g_count = self.sequence.count('G') + self.sequence.count('g')
        c_count = self.sequence.count('C') + self.sequence.count('c')

        gc_count = g_count + c_count

        return gc_count/total * 100

    def base_counts(self) -> dict[str, int]:
        """Return the count of each nucleotide symbol in the sequence,
            normalized to uppercase."""

        if self.molecule_type == MoleculeType.PROTEIN:
            raise ValueError('Cannot count bases in proteins')
        
        sequence = self.sequence.upper()
        counts = {}
        for letter in sequence:
            if letter in counts:
                counts[letter] += 1
            else:
                counts[letter] = 1
        return counts

    def base_frequencies(self) -> dict[str, float]:
        """Return the frequency of each symbol in the nucleotide sequence as a percentage.

            Raises:
                ValueError: If the sequence is empty or is a protein.
        """
        if self.molecule_type == MoleculeType.PROTEIN:
            raise ValueError('Cannot count bases in proteins')

        counts = self.base_counts()
        total = len(self.sequence)

        if total == 0:
            raise ValueError('Cannot calculate base frequencies of an empty sequence')

        frequencies = {}
        for base in counts:
            frequencies[base] = counts[base]/total * 100
        return frequencies

    def reverse_complement(self) -> str:
        """Return the reverse complement of the DNA or RNA sequence.

            Raises:
                ValueError: If the sequence contains an unsupported DNA base or is a protein.
        """
        if self.molecule_type == MoleculeType.PROTEIN:
            raise ValueError('Wrong Molecule')
        elif self.molecule_type == MoleculeType.DNA:
            complements = DNA_COMPLEMENTS
        else:
            complements = RNA_COMPLEMENTS
        reverse = self.sequence[::-1]
        rev_compl = ''
        for letter in reverse:
            if letter.upper() in complements:
                rev_compl += complements[letter.upper()]
            else:
                raise ValueError("Invalid base.")

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
            ValueError: If ``strand`` is not ``"coding"`` or ``"template"``
            or the molecule type is RNA or protein.
        """
        if self.molecule_type == MoleculeType.PROTEIN:
            raise ValueError('Proteins cannot do transcription')
        elif self.molecule_type == MoleculeType.RNA:
            raise ValueError('RNA cannot do transcription')

        if strand == 'coding':
            sequence = self.sequence.upper()
        elif strand == 'template':
            sequence = self.reverse_complement()
        else:
            raise ValueError('Invalid Strand')

        return sequence.replace('T', 'U')


    def translate(self, frame: int | None = None) -> str:     
        """
        Translate the sequence from the first start codon of DNA or RNA to the first stop codon.

        When ``frame`` is ``None``, the first start codon in the sequence is used.
        When ``frame`` is 1, 2, or 3, the first start codon in the selected
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
            ValueError: If ``frame`` is not ``None``, 1, 2, or 3.
                        If it is a protein.
         """
        if self.molecule_type == MoleculeType.PROTEIN:
            raise ValueError('Cannot translate Proteins')

        elif self.molecule_type == MoleculeType.RNA:
            sequence = self.sequence.replace('U','T')

        else:
            sequence = self.sequence.upper()

        if frame is None:
            start = sequence.find('ATG')

        else:
            if frame > 3 or frame < 1:
                raise ValueError('Invalid Frame')
            start = -1
            offset = frame - 1
            for i in range(offset,len(sequence)-2,3):
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
        """Return open reading frames found in the selected DNA or RNA strand.

        Searches for complete ORFs beginning with ``ATG`` and ending at a
        compatible stop codon. When `frame` is `None`, all reading frames are searched.
        When `frame` is 1, 2, or 3, the search is restricted to that reading frame.

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
                or ``"both"``, or if ``frame`` is not ``None``, 1, 2, or 3.
        """

        if self.molecule_type == MoleculeType.PROTEIN:
            raise ValueError('Invalid Molecule Type: Protein')

        results = []
        if strand == 'both':
            return self.find_orfs('forward',frame=frame) + self.find_orfs('reverse',frame=frame)
        if strand == "forward":
            sequence = self.sequence.upper()
            if self.molecule_type == MoleculeType.RNA:
                sequence = sequence.replace('U','T')
        elif strand == "reverse":
            sequence = self.reverse_complement()
            if self.molecule_type == MoleculeType.RNA:
                sequence = sequence.replace('U','T')
        else:
            raise ValueError('Invalid Strand')



        if frame is None:
            positions = range(len(sequence))
        elif 1 <= frame <= 3:
            offset = frame -1 
            positions = range(offset,len(sequence),3)
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
                        if self.molecule_type == MoleculeType.RNA:
                            orf = orf.replace('T','U')
                        results.append(orf)
                        break
                    else:
                        orf += triplet

        return results

    def kmers(self, k:int):
        """ Return all k-mers of length `k` from the sequence.

            Raises:
                ValueError: If invalid k number ( k <= 0).           
        """
        if k > 0:
            results = []
            kmer = ''
            for i in range(len(self.sequence)):
                kmer += self.sequence[i:i+k]
                if len(kmer) == k:
                    results.append(kmer)
                    kmer = ''
        else:
            raise ValueError('Invalid k number')
            
        return results

    def kmer_counts(self, k:int) -> dict:
        """Return all k-mers of length `k` from the sequence.
            Raises:
                ValueError: for k <= 0"""

        results = self.kmers(k)

        kmers = {}
        for kmer in results:
            freq = results.count(kmer)
            if kmer not in kmers:
                kmers[kmer] = freq

        return kmers

    def amino_acid_counts(self):
        """Return the count of each amino acid residue in the sequence, normalized to uppercase, restricted to proteins.
        
        Raises:
            ValueError: if molecule type is DNA or RNA.
        """

        if self.molecule_type == MoleculeType.DNA:
            raise ValueError('Cannot count aminoacids in DNA')
        elif self.molecule_type == MoleculeType.RNA:
            raise ValueError('Cannot count aminoacids in RNA')
        
        sequence = self.sequence.upper()
        counts = {}
        for letter in sequence:
            if letter in counts:
                counts[letter] += 1
            else:
                counts[letter] = 1

        return counts

    def amino_acid_frequencies(self):
        """Return the frequency of each amino acid residue in the sequence as a percentage, restricted to proteins.

            Raises:
                ValueError: If the sequence is empty or if molecule type is DNA or RNA.
        """
        counts = self.amino_acid_counts()
        total = len(self.sequence)

        if total == 0:
            raise ValueError('Cannot calculate amino acid frequencies of an empty sequence')

        frequencies = {}
        for aa in counts:
            frequencies[aa] = counts[aa]/total * 100
        return frequencies

    def filter_by_length(sequences, min_length = None, max_length = None):

        if min_length is not None and max_length is not None:
            if min_length > max_length:
                raise ValueError('Minimum length must not be higher than maximum length')

        elif min_length is not None:
            if min_length < 0:
                raise ValueError('Minimum length cannot be lower than 0')

        elif max_length is not None:
            if max_length < 0:
                raise ValueError('Maximum length cannot be lower than 0')

        else:
            raise ValueError('No length limits were provided')

        results = []
        for sequence in sequences:
            length = sequence.length()
            if min_length is not None and max_length is not None:
                if min_length <= length <= max_length:
                    results.append(sequence)
            elif min_length is not None:
                if length >= min_length:
                    results.append(sequence)
            elif max_length is not None:
                if length <= max_length:
                    results.append(sequence)

        return results

            


        


    