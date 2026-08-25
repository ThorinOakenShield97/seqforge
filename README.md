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
- Sequence length calculation
- Base counting and frequency analysis
- GC content calculation
- DNA reverse complement
- DNA → RNA transcription
- DNA → protein translation
- Translation from explicit reading frames
- IUPAC ambiguity-code support
- Motif searching
- Open reading frame (ORF) detection
- Forward, reverse, and both-strand ORF searches
- ORF searches by reading frame
- Command-line interface
- Python API
- Literal sequence and FASTA file input from the CLI
- CLI support for FASTA multi-record files

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

9
33.333333333333336
CTATTTCAT
AUGAAAUAG
MK

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

## Translation

By default, translation searches for the first start codon in the sequence:

```python
sequence.translate()
```

A specific reading frame can be selected with 0, 1, or 2:

```python
sequence.translate(frame=0)
sequence.translate(frame=1)
sequence.translate(frame=2)
```

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
sequence.find_orfs(frame=0)
sequence.find_orfs(strand="reverse", frame=1)
sequence.find_orfs(strand="both", frame=2)
```

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

## 💻 Command-line interface

SeqForge provides a command-line interface:

```bash
seqforge --help
```

Available commands:
    gc
    orf
    stats
    translate
    version

## Input

SeqForge commands accept either:

- a DNA sequence provided directly on the command line; or
- a path to an existing FASTA file.

For example:

```bash
seqforge stats ATGC
```
or

```bash
seqforge stats sequence.fasta
```

FASTA files may contain multiple records. Each record is processed independently and its FASTA identifier is preserved in the output.

## GC content

```bash
seqforge gc ATGC
```
GC content: 50.0%

A FASTA file can be used directly:

```bash
seqforge gc sequence.fasta
```

## Sequence statistics

```bash
seqforge stats ATGC
```

Example output:
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

For FASTA input, the record identifier is shown before its statistics.

## Translation

```bash
seqforge translate ATGAAATAG
```
Protein: MK

A reading frame can be selected with --frame:

```bash
seqforge translate AATGAAATAG --frame 1
```

FASTA files are also supported:

```bash
seqforge translate gene.fasta
```
Each FASTA record is translated independently and its identifier is preserved in the output.
A reading frame can be selected with --frame:

```bash
seqforge translate gene.fasta --frame 1
```

## ORF detection

```bash
seqforge orf ATGAAATAG
```

A reading frame can be selected with --frame:

```bash
seqforge orf AATGAAATAG --frame 1
```

ORFs can also be searched on different strands:

```bash
seqforge orf sequence.fasta --strand forward
seqforge orf sequence.fasta --strand reverse
seqforge orf sequence.fasta --strand both
```

When both strands are requested, the output is grouped explicitly:
    Forward:
    ATGTGA
    Reverse:
    ATGTAG

A frame can be combined with strand selection:

```bash
seqforge orf sequence.fasta --strand both --frame 1
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

### 0.2.0

The 0.2.0 release line focuses on expanding the command-line experience and improving input handling:

- CLI support for sequence and FASTA input
- FASTA multi-record processing
- Reading-frame selection for translation and ORF detection
- Strand-aware ORF analysis
- Consistent CLI input and error handling
- Sequence statistics CLI

### Future releases

Planned areas of development include:

- RNA sequence support
- DNA transcription from forward and reverse strands
- Protein sequence support
- FASTQ parsing and quality analysis
- Sequence filtering operations
- Additional sequence-analysis tools
- Additional biological sequence formats
- Expanded documentation
- Performance improvements
- Further API and CLI refinement

## 🤝 Contributing

Contributions, ideas, and bug reports are welcome.

Before submitting changes, please make sure the test suite passes:

```bash
uv run pytest
```

## 📄 License

SeqForge is licensed under the MIT License.