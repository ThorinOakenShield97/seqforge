from seqforge.models.sequence import Sequence

def pairwise_distances(sequences: list[Sequence]) -> list[tuple[str,str,int]]:
    """
    Calculate pairwise Hamming distances between sequences.

    Args:
        sequences: Sequences to compare pairwise.

    Returns:
        A list of tuples containing the IDs of each sequence pair and
        their Hamming distance.

    Raises:
        ValueError: If any pair of sequences has incompatible molecule types
        or different lengths.
    """

    results = []
    for i, seq_1 in enumerate(sequences):
        for seq_2 in sequences[i+1:]:
            distance = seq_1.distance(seq_2)

            result = (seq_1.id,seq_2.id, distance)
            results.append(result)

    return results