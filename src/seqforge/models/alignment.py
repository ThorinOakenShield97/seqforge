def global_alignment(seq1, seq2, match=1, mismatch=-1, gap=-1):

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

    return aligned_seq1, aligned_seq2