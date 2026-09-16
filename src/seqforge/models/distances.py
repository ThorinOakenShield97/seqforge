def pairwise_distances(sequences):

    results = []
    for i, seq_1 in enumerate(sequences):
        for seq_2 in sequences[i+1:]:
            distance = seq_1.distance(seq_2)

            result = (seq_1.id,seq_2.id, distance)
            results.append(result)

    return results