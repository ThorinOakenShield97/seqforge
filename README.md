# SeqForge

> A modern Python toolkit for biological sequence analysis.

Current release: 0.6.0

SeqForge is an open-source Python toolkit for common biological sequence analysis tasks.

It provides both a **Python API** for working with biological sequences and a **command-line interface (CLI)** for interacting with the project from the terminal.

The project focuses on clear APIs, predictable behaviour, automated testing, and modern Python development practices.

---

## ✨ Features

SeqForge currently provides:

- FASTA parsing
- FASTA multi-record processing
- FASTQ parsing
- FASTQ multi-record processing
- FASTQ quality analysis
- Sequence length calculation
- Base counting and frequency analysis for DNA and RNA
- GC content calculation for DNA and RNA
- GC content analysis by sliding windows
- DNA and RNA reverse complement
- DNA → RNA transcription
- RNA → complementary DNA reverse transcription
- Coding and template strand transcription
- DNA/RNA → protein translation
- Translation from explicit reading frames
- Biological reading frames 1, 2, and 3
- IUPAC ambiguity-code support
- Motif searching with overlapping matches
- 1-based motif positions
- Open reading frame (ORF) detection
- Forward, reverse, and both-strand ORF searches
- ORF searches by reading frame
- K-mer generation
- K-mer counting
- K-mer frequency analysis
- K-mer searching
- Sequence filtering by length
- Sequence filtering by motif
- FASTQ read filtering by minimum mean quality
- Hamming sequence distances
- Pairwise distance analysis
- Global sequence alignment
- Progressive multiple sequence alignment
- Command-line interface
- Python API
- Literal sequence, FASTA, and FASTQ input from the CLI
- DNA, RNA, and protein molecule-type support
- CLI support for FASTA and FASTQ multi-record files

---

## 🚀 Installation

SeqForge requires Python 3.12 or later.

### From source

Clone the repository and install it with `uv`:

```bash
uv sync
```

You can then run SeqForge with:

```bash
uv run seqforge --help
```

## 📖 Python API

The main API is provided by the `Sequence` class.

```python
from seqforge.models.sequence import Sequence

sequence = Sequence(
    id="example",
    sequence="ATGAAATAG",
)

print(sequence.length())
print(sequence.gc_content())
print(sequence.reverse_complement())
print(sequence.transcribe())
print(sequence.translate())
```

Output:

```text
9
33.333333333333336
CTATTTCAT
AUGAAAUAG
MK
```

## 🧬 Molecule types

SeqForge supports three explicit molecule types:

```python
from seqforge.models.molecule_type import MoleculeType
from seqforge.models.sequence import Sequence

dna = Sequence(
    id="dna",
    sequence="ATGC",
    molecule_type=MoleculeType.DNA,
)

rna = Sequence(
    id="rna",
    sequence="AUGC",
    molecule_type=MoleculeType.RNA,
)

protein = Sequence(
    id="protein",
    sequence="MKWVTF",
    molecule_type=MoleculeType.PROTEIN,
)
```

DNA is the default molecule type for backwards compatibility.

String values are also accepted case-insensitively:

```python
Sequence(id="rna", sequence="AUGC", molecule_type="RNA")
Sequence(id="protein", sequence="MKWVTF", molecule_type="protein")
```

Invalid molecule types raise `ValueError`.

## Sequence analysis

```python
sequence.base_counts()
```

```python
sequence.base_frequencies()
```

### Motif searching

Motifs are searched case-insensitively and overlapping matches are supported.

```python
sequence.find_motif("ATG")
```

Motif positions are reported using 1-based coordinates.

For example:

```python
Sequence(
    id="example",
    sequence="AAAA",
).find_motif("AA")
```

returns:

```text
[1, 2, 3]
```

A motif longer than the sequence returns an empty list.

## Filtering

SeqForge provides filtering operations for collections of sequences and FASTQ reads.

### Filter by length

```python
from seqforge.models.filters import filter_by_length

filtered = filter_by_length(
    sequences,
    min_length=100,
    max_length=500,
)
```

Minimum and maximum lengths are inclusive. Either limit may be omitted.

### Filter by motif

```python
from seqforge.models.filters import filter_by_motif

filtered = filter_by_motif(sequences, motif="ATG")
```

Only sequences containing the motif are retained.

### Filter FASTQ reads by quality

```python
from seqforge.models.filters import filter_by_quality

filtered = filter_by_quality(
    reads,
    min_quality=30,
)
```

Filtering is based on the read's mean Phred quality, using the quality metrics already provided by `FastqRead`.

## 🪟 GC content by sliding windows

GC content can be calculated over sliding windows:

```python
from seqforge.models.gc_window import gc_content_windows

windows = gc_content_windows(sequence, window_size=100)
```

Each result contains:

```text
(start, end, gc_content)
```

Positions are reported using 1-based coordinates.

For example:

```text
(1, 100, 48.0)
(2, 101, 49.0)
(3, 102, 50.0)
```

DNA and RNA sequences are supported. Protein sequences are not valid inputs for GC content analysis.

## Transcription

DNA sequences can be transcribed using coding or template strand semantics. RNA and protein sequences cannot be transcribed.

The default behaviour uses the coding strand:

```python
sequence.transcribe()
```

The coding strand can be selected explicitly:

```python
sequence.transcribe(strand="coding")
```

The template strand can also be transcribed:

```python
sequence.transcribe(strand="template")
```

The template strand is handled by obtaining its reverse complement before transcription.

## Reverse transcription

RNA sequences can be reverse transcribed into complementary DNA:

```python
sequence = Sequence(
    id="rna",
    sequence="AUGCCU",
    molecule_type="RNA",
)

result = sequence.reverse_transcribe()
print(result.sequence)
```

Output:

```text
TACGGA
```

The returned value is a new DNA `Sequence` preserving the original sequence identifier.

## Translation

DNA and RNA sequences can be translated. For RNA, SeqForge normalizes U-based codons internally.

By default, translation searches for the first start codon in the sequence:

```python
sequence.translate()
```

A specific reading frame can be selected with 1, 2, or 3.

```python
sequence.translate(frame=1)
sequence.translate(frame=2)
sequence.translate(frame=3)
```

Frame numbering is user-facing and follows biological convention:

- frame 1 → offset 0
- frame 2 → offset 1
- frame 3 → offset 2

## ORF detection

```python
sequence.find_orfs()
```

ORFs can be searched on the forward strand, reverse complement, or both:

```python
sequence.find_orfs(strand="forward")
sequence.find_orfs(strand="reverse")
sequence.find_orfs(strand="both")
```

A specific reading frame can also be selected:

```python
sequence.find_orfs(frame=1)
sequence.find_orfs(strand="reverse", frame=2)
sequence.find_orfs(strand="both", frame=3)
```

The six possible reading frames can therefore be represented as:

- +1 → forward, frame 1
- +2 → forward, frame 2
- +3 → forward, frame 3
- -1 → reverse, frame 1
- -2 → reverse, frame 2
- -3 → reverse, frame 3

## 🧬 IUPAC ambiguity codes

SeqForge supports IUPAC ambiguity codes when working with biological sequences.

For example:

- GCN

represents:

- GCA
- GCC
- GCG
- GCT

The `expand_iupac()` utility can be used directly:

```python
from seqforge.models.sequence import expand_iupac

expand_iupac("GCN")
```

IUPAC ambiguity is also supported by translation and ORF detection.

When all possible expansions of an ambiguous codon produce the same amino acid, that amino acid is used.

When the possible results differ, translation uses `X`.

## 🧬 K-mers

SeqForge supports k-mer generation, counting, frequency analysis, and searching.

K-mers are extracted using a sliding window:

```python
sequence.kmers(k=2)
```

For example:

```python
Sequence(
    id="example",
    sequence="ATGC",
).kmers(k=2)
```

returns:

```python
["AT", "TG", "GC"]
```

K-mer counts can be calculated with:

```python
sequence.kmer_counts(k=2)
```

For example:

```python
Sequence(
    id="example",
    sequence="ATAT",
).kmer_counts(k=2)
```

returns:

```python
{
    "AT": 2,
    "TA": 1,
}
```

K-mer frequencies can be calculated with:

```python
sequence.kmer_frequencies(k=3)
```

For example:

```python
Sequence(
    id="example",
    sequence="ATGATG",
).kmer_frequencies(k=3)
```

returns:

```python
{
    "ATG": 0.5,
    "TGA": 0.25,
    "GAT": 0.25,
}
```

Specific k-mers can be searched in a sequence with 1-based positions:

```python
sequence.find_kmers(["AT", "GC"])
```

The search is case-insensitive and overlapping matches are supported.

A non-positive k raises `ValueError`.

If k is greater than the sequence length, no k-mers are produced.

## 🧬 Sequence distances

SeqForge provides Hamming distance between two sequences:

```python
sequence_1 = Sequence(id="seq1", sequence="ATGC")
sequence_2 = Sequence(id="seq2", sequence="ATCC")

sequence_1.distance(sequence_2)
```

returns:

```text
1
```

The comparison is case-insensitive. Sequences must have equal lengths and the same molecule type.

Pairwise distances can be calculated for a collection of sequences:

```python
from seqforge.models.distances import pairwise_distances

pairwise_distances(sequences)
```

The result contains one tuple for each unique pair:

```text
[
    ("seq1", "seq2", 1),
    ("seq1", "seq3", 2),
    ("seq2", "seq3", 1),
]
```

## 🧬 Sequence alignment

SeqForge provides global alignment between two sequences:

```python
from seqforge.models.alignment import global_alignment

aligned = global_alignment(
    Sequence(id="seq1", sequence="ATGC"),
    Sequence(id="seq2", sequence="ATC"),
)
```

The function returns the two aligned strings. Match, mismatch, and gap scores can be configured:

```python
global_alignment(
    seq1,
    seq2,
    match=2,
    mismatch=-1,
    gap=-2,
)
```

Multiple sequence alignment uses a progressive strategy anchored on the first sequence:

```python
from seqforge.models.alignment import multiple_alignment

aligned = multiple_alignment(sequences)
```

The result is a list of aligned sequence strings with gaps represented by `-`.

---

## 🧪 FASTQ

SeqForge supports FASTQ input for read-oriented analysis.

A FASTQ record contains:

- a read identifier
- a nucleotide sequence
- a + separator
- a per-base quality string

SeqForge currently uses Phred+33 quality encoding.

For example:

```text
@read1
ATGC
+
IIII
```

represents a read with:

```text
ID:       read1
Sequence: ATGC
Quality:  IIII
```

FASTQ files may contain multiple records, which are processed independently. Molecule-aware analysis can apply DNA, RNA, or protein semantics to the sequence associated with each read.

`FastqRead` provides quality-related analysis:

```python
read.quality_scores()
read.mean_quality()
read.min_quality()
read.max_quality()
```

## 💻 Command-line interface

SeqForge provides a command-line interface:

```bash
seqforge --help
```

Available commands:

```text
align
distance
distances
filter
gc
kmer
motif
orf
reverse-transcribe
stats
transcribe
translate
version
```

## Input

SeqForge accepts:

- a biological sequence provided directly on the command line
- a path to an existing FASTA file
- a path to an existing FASTQ file

Not every command accepts every input format. FASTQ support is focused on commands that can operate on the nucleotide sequence of each read, together with dedicated FASTQ quality analysis where applicable.

Commands that support molecule-aware analysis accept `--molecule-type` with `DNA`, `RNA`, or `PROTEIN` (case-insensitive).

For example:

```bash
seqforge stats ATGC
```

```bash
seqforge stats sequence.fasta
```

```bash
seqforge stats sequence.fastq
```

FASTA and FASTQ files may contain multiple records. Each record is processed independently and its identifier is preserved in the output.

## Filtering from the CLI

The `filter` command supports length, motif, and FASTQ quality criteria.

Filter a FASTA file by minimum length:

```bash
seqforge filter sequences.fasta --min-length 100
```

Filter by maximum length:

```bash
seqforge filter sequences.fasta --max-length 500
```

Filter by motif:

```bash
seqforge filter sequences.fasta --motif ATG
```

Multiple criteria can be combined. Criteria are applied together:

```bash
seqforge filter sequences.fasta --min-length 100 --max-length 500 --motif ATG
```

FASTQ reads can be filtered by minimum mean quality:

```bash
seqforge filter reads.fastq --min-quality 30
```

The FASTQ records that pass filtering keep their identifier, sequence, separator, and quality string intact.

## GC content

```bash
seqforge gc ATGC
```

Output:

```text
GC content: 50.0%
```

A FASTA file can be used directly:

```bash
seqforge gc sequence.fasta
```

FASTQ input:

```bash
seqforge gc sequence.fastq
```

GC content is supported for DNA and RNA sequences. Protein sequences are not valid inputs for GC content analysis.

For FASTQ input, GC content is calculated from each read sequence and quality scores do not affect the calculation.

## GC content by window

Use `--window-size` to calculate GC content over sliding windows:

```bash
seqforge gc sequence.fasta --window-size 100
```

Output:

```text
>seq1
positions 1 to 100: 48.0%
positions 2 to 101: 49.0%
positions 3 to 102: 50.0%
```

Window positions are reported using 1-based coordinates. The same windowed analysis can be applied to FASTQ reads:

```bash
seqforge gc reads.fastq --window-size 100
```

## Protein statistics

When the input molecule type is protein, `stats` reports sequence length, amino acid counts, and amino acid frequencies instead of nucleotide statistics.

```bash
seqforge stats MKWVTFM --molecule-type PROTEIN
```

## Sequence statistics

```bash
seqforge stats ATGC
```

Example output:

```text
Length: 4
GC content: 50.0%
Base counts:
A: 1
C: 1
G: 1
T: 1
Base frequencies:
A: 25.0%
C: 25.0%
G: 25.0%
T: 25.0%
```

For FASTA input, the record identifier is shown before its statistics.

## FASTQ statistics

`stats` provides additional read-oriented metrics for FASTQ files.

Example:

```bash
seqforge stats reads.fastq
```

Output:

```text
@read1
Length: 150
GC content: 48.0%
Mean quality: 32.4
Min quality: 18
Max quality: 40
@read2
Length: 150
GC content: 51.3%
Mean quality: 31.8
Min quality: 20
Max quality: 40
Reads: 2
Mean read length: 150.0
Overall mean quality: 32.1
Mean GC content: 49.65%
Min read length: 150
Max read length: 150
```

## Motif command

Search for a motif in a literal sequence:

```bash
seqforge motif ATGCATGC ATGC
```

Output:

```text
positions: [1, 5]
```

Motif positions use 1-based coordinates and overlapping matches are supported.

FASTA and FASTQ input are also supported:

```bash
seqforge motif sequences.fasta ATG
```

```bash
seqforge motif reads.fastq ATG
```

Each record is processed independently.

## Transcription

```bash
seqforge transcribe ATGC
```

Output:

```text
AUGC
```

The strand can be selected explicitly:

```bash
seqforge transcribe GCAT --strand template
```

Both coding and template strands can be requested:

```bash
seqforge transcribe ATGC --strand both
```

Output:

```text
Coding:
AUGC
Template:
GCAU
```

FASTA and FASTQ input are supported:

```bash
seqforge transcribe gene.fasta --strand template
```

```bash
seqforge transcribe reads.fastq --strand coding
```

Each record is processed independently.

## Reverse transcription

Reverse transcribe an RNA sequence into complementary DNA:

```bash
seqforge reverse-transcribe AUGCCU
```

Output:

```text
TACGGA
```

FASTA and FASTQ input are also supported:

```bash
seqforge reverse-transcribe reads.fasta
```

```bash
seqforge reverse-transcribe reads.fastq
```

## Translation

```bash
seqforge translate ATGAAATAG
```

```text
Protein: MK
```

A reading frame can be selected with `--frame`:

```bash
seqforge translate AATGAAATAG --frame 2
```

FASTA and FASTQ files are also supported:

```bash
seqforge translate gene.fasta
```

```bash
seqforge translate reads.fastq --molecule-type RNA
```

Each record is translated independently and its identifier is preserved in the output where applicable.

Reading frames are numbered 1, 2, and 3.

## ORF detection

```bash
seqforge orf ATGAAATAG
```

A reading frame can be selected with `--frame`:

```bash
seqforge orf AATGAAATAG --frame 2
```

ORFs can also be searched on different strands:

```bash
seqforge orf sequence.fasta --strand forward
seqforge orf sequence.fasta --strand reverse
seqforge orf sequence.fasta --strand both
```

When both strands are requested, the output is grouped explicitly:

```text
Forward:
ATGTGA
Reverse:
ATGTAG
```

A frame can be combined with strand selection:

```bash
seqforge orf sequence.fasta --strand both --frame 1
```

FASTQ input is also supported:

```bash
seqforge orf reads.fastq --strand reverse
```

## K-mer command

K-mers can be extracted from DNA, RNA, or protein sequences provided as literals, FASTA files, or FASTQ files.

```bash
seqforge kmer ATGC --k 2
```

Output:

```text
AT
TG
GC
```

K-mer counts can be requested with `--counts`:

```bash
seqforge kmer ATAT --k 2 --counts
```

Output:

```text
AT: 2
TA: 1
```

K-mer frequencies can be requested with `--frequencies`:

```bash
seqforge kmer ATGATG --k 3 --frequencies
```

K-mers can also be searched with `--find`:

```bash
seqforge kmer ATGATG --find ATG --find TGA
```

The `--counts`, `--frequencies`, and `--find` modes are mutually exclusive.

FASTA input is supported:

```bash
seqforge kmer sequence.fasta --k 3
```

FASTQ input is also supported:

```bash
seqforge kmer reads.fastq --k 3
```

For FASTQ input, k-mer analysis is performed on the nucleotide sequence of each read; quality information is not used by the basic k-mer command.

## Sequence distance

Calculate the Hamming distance between two sequences:

```bash
seqforge distance ATGC ATCC
```

Output:

```text
distance: 1
```

The same command supports RNA and protein sequences by using `--molecule-type`:

```bash
seqforge distance AUGC AUCC --molecule-type RNA
```

```bash
seqforge distance MKWVTF MKAVTF --molecule-type PROTEIN
```

The two sequences must have the same length and molecule type.

## Pairwise distances

The `distances` command calculates all unique pairwise distances from a FASTA or FASTQ collection:

```bash
seqforge distances sequences.fasta
```

Example output:

```text
seq1, seq2: 1
seq1, seq3: 2
seq2, seq3: 1
```

FASTQ input is supported:

```bash
seqforge distances reads.fastq
```

## Sequence alignment

Align two or more sequences from a literal input file:

```bash
seqforge align sequences.fasta
```

For a multiple FASTA file, aligned sequences are printed one per line:

```text
ATGC
AT-C
ATGC
```

FASTQ input is also supported:

```bash
seqforge align reads.fastq
```

RNA and protein alignment can be requested explicitly:

```bash
seqforge align sequences.fasta --molecule-type RNA
```

```bash
seqforge align proteins.fasta --molecule-type PROTEIN
```

## 🧪 Development

SeqForge uses uv for project and dependency management.

Run the test suite with:

```bash
uv run pytest
```

Build the package with:

```bash
uv build
```

The project maintains an automated test suite covering core sequence-analysis functionality, FASTA parsing, FASTQ parsing, CLI behaviour, filtering, windowed GC analysis, motif analysis, k-mer analysis, sequence distances, sequence alignment, and package functionality.

GitHub Actions run the test suite and verify that the package can be built.

## 🎯 Project philosophy

SeqForge follows one simple idea:

Do one thing and do it well.

The project aims to provide:

- clear and predictable APIs
- reproducible results
- modern Python practices
- strong automated testing
- lightweight and modular components
- useful bioinformatics functionality without unnecessary complexity
- an intuitive command-line experience

## 🗺️ Roadmap


### 0.5.0

The 0.5.x release line focused on sequence filtering and windowed analysis:

- Sequence filtering by length
- Sequence filtering by motif
- FASTQ filtering by minimum mean quality
- CLI filtering with combined criteria
- GC content analysis by sliding windows
- Windowed GC analysis for FASTA and FASTQ input
- 1-based window coordinates in user-facing results
- Expanded CLI and API documentation

### 0.6.0

The 0.6.x release line focuses on deeper sequence analysis:

- Motif analysis with 1-based positions and overlapping matches
- RNA → complementary DNA reverse transcription
- K-mer frequency analysis
- K-mer searching
- Hamming sequence distances
- Pairwise distance analysis
- Global sequence alignment
- Progressive multiple sequence alignment
- Expanded FASTQ support across sequence-analysis commands
- Expanded CLI and API documentation

### Future releases

Planned areas of development include:

- Additional biological sequence formats
- GFF/GFF3 annotation support
- Region-based sequence analysis
- Performance improvements
- Streaming support for large files
- Further API and CLI refinement

## 🤝 Contributing

Contributions, ideas, and bug reports are welcome.

Before submitting changes, please make sure the test suite passes:

```bash
uv run pytest
```

## 📄 License

SeqForge is licensed under the MIT License.
