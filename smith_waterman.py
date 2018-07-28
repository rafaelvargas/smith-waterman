# SEQUENCE ALIGNMENT VALUES TO DETERMINE SCORE
MATCH = 1
MISMATCH = -1
GAP = -2


# Initialize score and traceback matrices
def matrix_initialization(seq1, seq2):
    nbr_lines = len(seq1) + 1
    nbr_columns = len(seq2) + 1
    score_matrix = [[0] * nbr_columns for i in range(nbr_lines)]
    traceback_matrix = [[(0, 0)] * nbr_columns for i in range(nbr_lines)]

    # Traceback matrix initialization
    for i in range(1, nbr_columns):
        traceback_matrix[0][i] = (0, i - 1)

    for i in range(1, nbr_lines):
        traceback_matrix[i][0] = (i - 1, 0)

    return score_matrix, traceback_matrix


def scoring(seq1, seq2, score_matrix, traceback_matrix):
    partial = 0  # diagonal comparing result

    for i in range(1, len(seq1) + 1):
        for j in range(1, len(seq2) + 1):

            if seq1[i - 1] == seq2[j - 1]:
                partial = MATCH
            else:
                partial = MISMATCH

            up = score_matrix[i - 1][j] + GAP
            left = score_matrix[i][j - 1] + GAP
            diag = score_matrix[i - 1][j - 1] + partial
            score_matrix[i][j] = max(up, left, diag, 0)

            if score_matrix[i][j] == up:
                traceback_matrix[i][j] = (i - 1, j)
            else:
                if score_matrix[i][j] == diag or score_matrix[i][j] == 0:
                    traceback_matrix[i][j] = (i - 1, j - 1)
                else:
                    traceback_matrix[i][j] = (i, j - 1)

    return score_matrix, traceback_matrix


# Find highest score in score matrix
def find_highest_score(score_matrix, nbr_lines, nbr_columns):
    highest_score = 0
    line, column = 0, 0
    for i in range(0, nbr_lines):
        for j in range(0, nbr_columns):
            if score_matrix[i][j] > highest_score:
                highest_score = score_matrix[i][j]
                line = i
                column = j

    return highest_score, (line, column)


# Also known as traceback step
def alignment(seq1, seq2, highest_score_pos, score_matrix, traceback_matrix):
    aligned_seq1 = ""
    aligned_seq2 = ""
    i = highest_score_pos[0]
    j = highest_score_pos[1]

    while score_matrix[i][j] != 0:
        if traceback_matrix[i][j][0] == i - 1 and traceback_matrix[i][j][1] == j - 1:
            aligned_seq1 += seq1[i - 1]
            aligned_seq2 += seq2[j - 1]
            i -= 1
            j -= 1
        else:
            if traceback_matrix[i][j][0] == i - 1:  # lines
                aligned_seq1 += seq1[i - 1]
                aligned_seq2 += "-"
                i -= 1
            else:  # columns
                aligned_seq1 += "-"
                aligned_seq2 += seq2[j - 1]
                j -= 1

    aligned_seq1 = aligned_seq1[::-1]
    aligned_seq2 = aligned_seq2[::-1]

    return aligned_seq1, aligned_seq2


def smith_waterman(seq1, seq2):
    score_matrix, traceback_matrix = matrix_initialization(seq1, seq2)
    score_matrix, traceback_matrix = scoring(seq1, seq2, score_matrix, traceback_matrix)
    highest_score, highest_score_pos = find_highest_score(
        score_matrix, len(seq1), len(seq2)
    )
    aligned_seq1, aligned_seq2 = alignment(
        seq1, seq2, highest_score_pos, score_matrix, traceback_matrix
    )

    return score_matrix, highest_score, aligned_seq1, aligned_seq2

