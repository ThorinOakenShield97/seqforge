from seqforge.models.sequence import Sequence

def global_alignment(seq1: Sequence, seq2: Sequence, match:int =1, mismatch:int =-1, gap:int =-1) -> tuple[str,str]:
    """Globally align two biological sequences.

    Args:
        seq1: First sequence to align.
        seq2: Second sequence to align.
        match: Score assigned to matching symbols.
        mismatch: Score assigned to mismatching symbols.
        gap: Score assigned to introducing a gap.

    Returns:
        A tuple containing the two aligned sequences.

    Raises:
        ValueError: If the sequences have different molecule types.
    """
    
    if seq1.molecule_type != seq2.molecule_type:
        raise ValueError('Sequences must have the same molecule type')

    def align(seq1,seq2,i,j,match, mismatch, gap, memo):
        if i == len(seq1.sequence):
            memo[(i, j)] = (len(seq2.sequence) - j) * gap
            return memo[(i, j)]

        if j == len(seq2.sequence):
            memo[(i, j)] = (len(seq1.sequence) - i) * gap
            return memo[(i, j)]

        if (i,j) in memo:
            return memo[(i,j)]

        if seq1.sequence[i] == seq2.sequence[j]:
            score = match

        else:
            score = mismatch

        diagonal_score = align(seq1,seq2,i+1,j+1,match, mismatch,gap,memo) + score
        gap_seq2 = align(seq1, seq2, i + 1, j, match, mismatch, gap, memo) + gap
        gap_seq1 = align(seq1, seq2, i, j+1, match, mismatch, gap, memo) + gap

        result = max(diagonal_score,gap_seq1,gap_seq2)

        memo[(i,j)] = result
        return memo[(i,j)]

    memo = dict()
    align(seq1, seq2, 0, 0, match, mismatch, gap, memo)

    i = 0
    j = 0

    aligned_seq1 = ''
    aligned_seq2 = ''

    while i < len(seq1.sequence) and j < len(seq2.sequence):
        if seq1.sequence[i] == seq2.sequence[j]:
            score = match
        else:
            score = mismatch

        if memo[(i, j)] ==  memo[(i + 1, j + 1)] + score:
            aligned_seq1 += seq1.sequence[i]
            aligned_seq2 += seq2.sequence[j]
            i += 1
            j += 1
        elif memo[(i, j)] == memo[(i + 1, j)] + gap:
            aligned_seq1 += seq1.sequence[i]
            aligned_seq2 += '-'
            i += 1
        elif memo[(i,j)] == memo[(i,j+1)] + gap:
            aligned_seq1 += '-'
            aligned_seq2 += seq2.sequence[j]
            j+= 1

    while i < len(seq1.sequence):
        aligned_seq1 += seq1.sequence[i]
        aligned_seq2 += '-'
        i += 1

    while j < len(seq2.sequence):
        aligned_seq1 += '-'
        aligned_seq2 += seq2.sequence[j]
        j += 1

    return aligned_seq1, aligned_seq2

def multiple_alignment(sequences: list[Sequence]) -> list[str]:
    """Align multiple biological sequences progressively.

    Args:
        sequences: Sequences to align. The first sequence is used as the
        reference for progressive alignment.

    Returns:
        A list containing the aligned sequences.

    Raises:
        ValueError: If no sequences are provided or the sequences have
        incompatible molecule types.
    """

    if not sequences:
        raise ValueError('No sequences in list to be aligned')

    if len(sequences) == 1:
        return [sequences[0].sequence]

    aligned = global_alignment(sequences[0],sequences[1])
    aligned = list(aligned)

    for seq in sequences[2:]:
        aligned_reference, aligned_new = global_alignment(sequences[0], seq)

        i = 0
        j = 0

        while i < len(aligned[0]) and j < len(aligned_reference):
            if aligned[0][i] == aligned_reference[j]:
                i += 1
                j += 1   

            elif aligned_reference[j] == '-':
                insert_pos = i
                for k in range(len(aligned)):
                    aligned[k] = (aligned[k][:insert_pos] + '-' + aligned[k][insert_pos:])
                i += 1
                j += 1 

            elif aligned[0][i] == '-':
                aligned_new = aligned_new[:j] + '-' + aligned_new[j:]
                i += 1

        while j < len(aligned_reference):
            if aligned_reference[j] == '-':
                for k in range(len(aligned)):
                    aligned[k] = aligned[k] + '-'
            j += 1

        aligned.append(aligned_new)   

    return aligned