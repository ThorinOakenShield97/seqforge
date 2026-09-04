from seqforge.models.sequence import Sequence

def gc_content_windows(sequence, window_size: int):
    """Calculate GC content over sliding windows.

    Args:
        sequence: DNA or RNA sequence to analyze.
        window_size: Size of each sliding window.

    Returns:
        A list of tuples containing the 1-based start position,
        1-based end position, and GC content of each window.

    Raises:
        ValueError: If window_size is not greater than 0.
    """

    if window_size is not None:
        if window_size <= 0:
            raise ValueError('Window must be greater than 0')

    results = []
    start = 1
    end = window_size
    length = sequence.length()

    while end <= length:
        curr = Sequence(id = sequence.id, sequence = sequence.sequence[start-1:end], molecule_type = sequence.molecule_type)
        gc = curr.gc_content()
        result = (start,end,gc)
        results.append(result)
        start += 1
        end += 1

    return results