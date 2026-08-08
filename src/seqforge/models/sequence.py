from dataclasses import dataclass
from itertools import product


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

def expand_iupac(codon):
    possibilities = [IUPAC_BASES[letter] for letter in codon]
    return [''.join(combination) for combination in product(*possibilities)]

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

    def base_counts(self):
        sequence = self.sequence.upper()
        counts = {}
        for letter in sequence:
            if letter in counts:
                counts[letter] += 1
            else:
                counts[letter] = 1
        return counts

    def base_frequencies(self):
        counts = self.base_counts()
        total = len(self.sequence)

        if total == 0:
            raise ValueError('Cannot calculate base frequencies of an empty sequence')

        frequencies = {}
        for base in counts:
            frequencies[base] = counts[base]/total * 100
        return frequencies

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


    def translate(self):
        sequence = self.sequence.upper()
        start = sequence.find('ATG')

        if start != -1:
            protein = ''
            for i in range(start,len(sequence)-2,3):
                triplet = sequence[i:i+3]
                if triplet in CODON_TABLE:
                    aa = CODON_TABLE[triplet]
                    if aa == '*':
                        break
                    protein += aa
                else:
                    possibilities = expand_iupac(triplet)
                    aminoacids = [CODON_TABLE[codon] for codon in possibilities]
                    if len(set(aminoacids)) == 1:
                        protein += aminoacids[0]
                    else:
                        protein += 'X'
        else:
            return ''
                
        return protein

    def find_motif(self,motif):
        positions = []
        sequence = self.sequence.upper()
        motif = motif.upper()
        if len(motif) == 0:
            raise ValueError('Empty Motif')
        for i in range(len(sequence)):
            if sequence[i:i+len(motif)] == motif:
                positions.append(i)

        return positions

    def find_orfs(self):
        results = []
        sequence = self.sequence.upper()
        for i in range(len(sequence)):
            triplet = sequence[i:i+3]
            if triplet == 'ATG':
                orf = ''
                orf += triplet
                for j in range(i+3,len(sequence),3):
                    triplet = sequence[j:j+3]
                    if triplet in CODON_TABLE:
                        if CODON_TABLE[triplet] == '*':
                            orf += triplet
                            results.append(orf)
                            break
                        orf += triplet
        return results
