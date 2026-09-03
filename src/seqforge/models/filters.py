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

def filter_by_motif(sequences, motif: str):

    results = []
    for sequence in sequences:
        result = sequence.find_motif(motif)
        if result:
            results.append(sequence)
    return results