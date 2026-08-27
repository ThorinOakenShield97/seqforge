# SeqForge

> A modern Python toolkit for biological sequence analysis.

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
- Base counting and frequency analysis
- GC content calculation
- DNA reverse complement
- DNA → RNA transcription
- Coding and template strand transcription
- DNA → protein translation
- Translation from explicit reading frames
- Biological reading frames 1, 2, and 3
- IUPAC ambiguity-code support
- Motif searching
- Open reading frame (ORF) detection
- Forward, reverse, and both-strand ORF searches
- ORF searches by reading frame
- K-mer generation
- K-mer counting
- Command-line interface
- Python API
- Literal sequence, FASTA, and FASTQ input from the CLI
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

The main API is provided by the Sequence class.

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

## Sequence analysis

```python
sequence.base_counts()
```
```python
sequence.base_frequencies()
```
```python
sequence.find_motif("ATG")
```

## Transcription

DNA sequences can be transcribed using coding or template strand semantics.

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

## Translation

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

The expand_iupac() utility can be used directly:

```python
from seqforge.models.sequence import expand_iupac

expand_iupac("GCN")
```

IUPAC ambiguity is also supported by translation and ORF detection.

When all possible expansions of an ambiguous codon produce the same amino acid, that amino acid is used.

When the possible results differ, translation uses 'X'.

## 🧬 K-mers

SeqForge supports k-mer generation and counting.

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
K-mer frequencies can be calculated with:

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

A non-positive k raises ValueError.
If k is greater than the sequence length, no k-mers are produced.

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
``

FASTQ files may contain multiple records, which are processed independently.

FastqRead provides quality-related analysis:

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
    gc
    kmer
    orf
    stats
    transcribe
    translate
    version
``
## Input

SeqForge accepts:

- a DNA sequence provided directly on the command line
- a path to an existing FASTA file
- a path to an existing FASTQ file

Not every command accepts every input format. FASTQ support is currently focused on read-oriented analysis such as statistics, GC content, and k-mer analysis.

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

## GC content

```bash
seqforge gc ATGC
```
Output:
GC content: 50.0%

A FASTA file can be used directly:

```bash
seqforge gc sequence.fasta
```

FASTQ input:

```bash
seqforge gc sequence.fastq
```

For FASTQ input, GC content is calculated from each read sequence and quality scores do not affect the calculation.

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

stats provides additional read-oriented metrics for FASTQ files.

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
``

## Transcription

```bash
seqforge transcribe ATGC
```
Output:
AUGC

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

FASTA input is supported:

```bash
seqforge transcribe gene.fasta --strand template
```
Each FASTA record is processed independently.

## Translation

```bash
seqforge translate ATGAAATAG
```
Protein: MK

A reading frame can be selected with --frame:

```bash
seqforge translate AATGAAATAG --frame 2
```

FASTA files are also supported:

```bash
seqforge translate gene.fasta
```
Each FASTA record is translated independently and its identifier is preserved in the output.
Reading frames are numbered 1, 2, and 3.

```bash
seqforge translate gene.fasta --frame 1
```

## ORF detection

```bash
seqforge orf ATGAAATAG
```

A reading frame can be selected with --frame:

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

## K-mer command

K-mers can be extracted from literal sequences, FASTA files, and FASTQ files.

```bash
seqforge kmer ATGC --k 2
```
Output:
```text
AT
TG
GC
```
K-mer frequencies can be requested with --counts:

```bash
seqforge kmer ATAT --k 2 --counts
```
Output:
```text
AT: 2
TA: 1
```

FASTA input is supported:

```bash
seqforge kmer sequence.fasta --k 3
```

FASTQ input is also supported:

```bash
seqforge kmer reads.fastq --k 3
```

For FASTQ input, k-mer analysis is performed on the nucleotide sequence of each read; quality information is not used by the basic k-mer command.

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

The project maintains an automated test suite covering core sequence-analysis functionality, FASTA parsing, CLI behaviour, and package functionality.

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

### 0.3.0

The 0.3.x release line focuses on expanding sequence analysis and introducing FASTQ workflows:

- Biological reading frames 1, 2, and 3
- Strand-aware ORF analysis
- Coding and template strand transcription
- FASTQ parsing
- FASTQ multi-record processing
- FASTQ quality analysis
- FASTQ statistics and read metrics
- K-mer generation and counting
- FASTA and FASTQ multi-record CLI processing
- Expanded CLI functionality

### Future releases

Planned areas of development include:

- RNA sequence support
- Protein sequence support
- Molecule-type aware sequence models
- Sequence filtering operations
- Quality-based FASTQ filtering
- GC content analysis by window
- Motif analysis improvements
- Additional k-mer analysis
- Simple sequence alignments
- Additional biological sequence formats
- GFF/GFF3 annotation support
- Region-based sequence analysis
- Performance improvements
- Streaming support for large files
- Expanded documentation
- Further API and CLI refinement

## 🤝 Contributing

Contributions, ideas, and bug reports are welcome.

Before submitting changes, please make sure the test suite passes:

```bash
uv run pytest
```

## 📄 License

SeqForge is licensed under the MIT License.
