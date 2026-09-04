def filter_by_length(sequences, min_length = None, max_length = None):
    """Filter sequences by their length.

    Args:
        sequences: Collection of sequences to filter.
        min_length: Minimum accepted sequence length, inclusive.
        max_length: Maximum accepted sequence length, inclusive.

    Returns:
        A list containing the sequences that satisfy the specified length limits.

    Raises:
        ValueError: If no length limit is provided, if a limit is negative,
            or if min_length is greater than max_length.
    """
        
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

def filter_by_motif(sequences, motif: str):

    """Filter sequences that contain a given motif.

    Args:
        sequences: Collection of sequences to filter.
        motif: Motif that must be present in each accepted sequence.

    Returns:
        A list containing the sequences that contain the specified motif.

    Raises:
        ValueError: If motif is empty.
    """

    results = []
    for sequence in sequences:
        result = sequence.find_motif(motif)
        if result:
            results.append(sequence)
    return results

def filter_by_quality(sequences, min_quality:int):
    """Filter FASTQ reads by their minimum mean quality.

    Args:
        sequences: Collection of FASTQ reads to filter.
        min_quality: Minimum accepted mean quality score, inclusive.

    Returns:
        A list containing the reads whose mean quality meets or exceeds
        the specified threshold.

    Raises:
        ValueError: If min_quality is None or negative.
    """

    if min_quality is not None:
        if min_quality < 0:
            raise ValueError('Minimum quality cannot be lower than 0')
    elif min_quality is None:
        raise ValueError('Minimum quality must be specified')

    results = []
    for sequence in sequences:
        quality = sequence.mean_quality()

        if quality >= min_quality:
            results.append(sequence)

    return results
    