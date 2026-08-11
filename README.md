# SeqForge

> A modern Python toolkit for biological sequence analysis.

SeqForge is an open-source Python toolkit for common biological sequence analysis tasks.

It provides both a **Python API** for working with biological sequences and a **command-line interface (CLI)** for interacting with the project from the terminal.

The project focuses on clear APIs, predictable behaviour, automated testing, and modern Python development practices.

---

## ✨ Features

SeqForge currently provides:

- FASTA parsing
- Sequence length calculation
- Base counting and frequency analysis
- GC content calculation
- DNA reverse complement
- DNA → RNA transcription
- DNA → protein translation
- IUPAC ambiguity-code support
- Motif searching
- Open reading frame (ORF) detection
- Forward, reverse, and both-strand ORF searches
- Command-line interface
- Python API

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

Example:
```text
SeqForge
Version: 0.1.0
Python: 3.14.6
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

The project maintains an automated test suite covering the core sequence-analysis functionality and CLI.

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

## 🗺️ Roadmap

### 0.1.x

The initial release focuses on core DNA sequence analysis:

- FASTA parsing
- Sequence analysis
- GC content
- Base statistics
- Reverse complement
- Transcription
- Translation
- IUPAC support
- Motif searching
- ORF detection
- CLI
- Package distribution

### Future releases

Planned areas of development include:

- FASTQ support
- Translation from explicit reading frames
- Additional sequence filtering operations
- File validation
- Expanded CLI functionality
- More comprehensive documentation
- Continuous integration
- Additional sequence-analysis tools

## 🤝 Contributing

Contributions, ideas, and bug reports are welcome.

Before submitting changes, please make sure the test suite passes:

```bash
uv run pytest
```

## 📄 License

SeqForge is licensed under the MIT License.