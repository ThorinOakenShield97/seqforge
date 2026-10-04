import typer

from seqforge.models.sequence import Sequence
from seqforge.commands.input import InputSource, resolve_input
from seqforge.exceptions import InvalidFastaError, InvalidFastqError

def transcribe(sequence:str, strand: str | None = None, molecule_type: str = 'dna') -> None:
    """Transcribe a DNA sequence into RNA.

    Args:
        sequence: Literal sequence or path to a FASTA or FASTQ file.
        strand: Strand to transcribe: ``"coding"``, ``"template"``,
            or ``"both"``. ``None`` uses the coding strand.
        molecule_type: Molecule type of the input sequence. Defaults to DNA.

    Raises:
        InvalidFastaError: If the FASTA input is invalid.
        InvalidFastqError: If the FASTQ input is invalid.
        FileNotFoundError: If the input file does not exist.
        ValueError: If the molecule type, strand, or sequence is invalid,
        or if the input is RNA or protein.
    """
    try:
        results = resolve_input(sequence, molecule_type)

        if results.source == InputSource.FASTQ:
            records = [Sequence(id=record.id, sequence=record.sequence, molecule_type=molecule_type) for record in results.records]
        else:
            records = results.sequences

        for seq in records:
            if results.source == InputSource.FASTA:
                print(f">{seq.id}")

            elif results.source == InputSource.FASTQ:
                print(f"@{seq.id}")

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
    except InvalidFastqError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code=1)
    except ValueError as e:
        typer.echo(f"Error: {e}", err=True)
        raise typer.Exit(code=1)
