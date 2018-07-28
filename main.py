from sys import argv

from file_io import read_sequences, write_matrix, write_alignments
from smith_waterman import smith_waterman


def main(argv):
    filename = argv[1]
    sequences = read_sequences(filename)

    print("Sequences in file:")
    for i, s in zip(range(1, len(sequences) + 1), sequences):
        print(str(i) + ": " + s[0])

    print("\nChoose the index (n,m) of the sequences that will be compared:")
    n = int(input("n: ")) - 1
    m = int(input("m: ")) - 1

    score_matrix, highest_score, aligned_seq1, aligned_seq2 = smith_waterman(
        sequences[n][1], sequences[m][1]
    )

    write_matrix(score_matrix, n, m)
    write_alignments(aligned_seq1, aligned_seq2, n, m, highest_score)


# Command line arguments: filename
if __name__ == "__main__":
    main(argv)

