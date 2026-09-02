import typer

from seqforge.models.sequence import Sequence
from seqforge.commands.input import InputSource, resolve_input
from seqforge.exceptions import InvalidFastaError

def transcribe(sequence:str, strand: str | None = None, molecule_type: str = 'dna') -> None:
    """Transcribe DNA sequences into RNA.

    Args:
        sequence: Literal sequence or path to a FASTA/FASTQ file.
        strand: Strand to transcribe: ``"coding"``, ``"template"``,
            or ``"both"``. ``None`` uses the coding strand.
        molecule_type: Molecule type of the input sequence. Defaults to DNA.

    Raises:
        ValueError: If the molecule type or strand is invalid, or if the input sequence is RNA or protein.
        InvalidFastaError: If the FASTA input is invalid.
        FileNotFoundError: If the input file does not exist.
    """
    try:
        results = resolve_input(sequence, molecule_type)
        for seq in results.sequences:
            if results.source == InputSource.FASTA:
                print(f">{seq.id}")

            if strand == 'both':
                coding = seq.transcribe(strand = 'coding')
                template = seq.transcribe(strand = 'template')
                print('Coding:')
                print(coding)
                print('Template:')        
                print(template)
            else:
                if strand is None:
                    rna = seq.transcribe()
                else:
                    rna = seq.transcribe(strand = strand)
                print(rna)
            
    except InvalidFastaError as i:
        typer.echo(f"Error: {i}", err=True)
        raise typer.Exit(code=1)
    except FileNotFoundError as f:
        typer.echo(f"Error: {f}", err=True)
        raise typer.Exit(code=1)
