#!/usr/bin/env python

def needleman_wunsch(seq_a: str, seq_b: str, match: int, mismatch: int, gap: int) -> tuple[tuple[str, str], int]:
    n = len(seq_a)
    m = len(seq_b)

    # Initialize the score matrix with zeros
    score_matrix = [[0] * (m + 1) for _ in range(n + 1)]

    # Initialize the traceback matrix
    traceback = [[0] * (m + 1) for _ in range(n + 1)]

    # Initialize the first row and column of the score matrix for gaps
    for i in range(1, n + 1):
        score_matrix[i][0] = score_matrix[i - 1][0] + gap
        traceback[i][0] = 1

    for j in range(1, m + 1):
        score_matrix[0][j] = score_matrix[0][j - 1] + gap
        traceback[0][j] = 2

    # Fill in the score matrix
    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match_mismatch = score_matrix[i - 1][j - 1] + (match if seq_a[i - 1] == seq_b[j - 1] else mismatch)
            insert_gap = score_matrix[i - 1][j] + gap
            delete_gap = score_matrix[i][j - 1] + gap
            score_matrix[i][j] = max(match_mismatch, insert_gap, delete_gap)

            if score_matrix[i][j] == match_mismatch:
                traceback[i][j] = 0  # Diagonal
            elif score_matrix[i][j] == insert_gap:
                traceback[i][j] = 1  # Up
            else:
                traceback[i][j] = 2  # Left

    # Traceback to find the alignment
    aligned_a, aligned_b = [], []
    i, j = n, m

    while i > 0 or j > 0:
        if traceback[i][j] == 0:  # Diagonal
            aligned_a.append(seq_a[i - 1])
            aligned_b.append(seq_b[j - 1])
            i -= 1
            j -= 1
        elif traceback[i][j] == 1:  # Up
            aligned_a.append(seq_a[i - 1])
            aligned_b.append('-')
            i -= 1
        else:  # Left
            aligned_a.append('-')
            aligned_b.append(seq_b[j - 1])
            j -= 1

    # Reverse the aligned sequences to get the correct order
    aligned_a = ''.join(aligned_a[::-1])
    aligned_b = ''.join(aligned_b[::-1])

    return (aligned_a, aligned_b), score_matrix[n][m]


