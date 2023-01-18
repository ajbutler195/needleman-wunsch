import sys

# store the score for an indel to be used.
indel = -1

# define a score function to determine the match score.
def S(A, B):
    if A == B:
        return 1
    else:
        return -1

def align(s1, s2):
    # initialize an array to be used for scoring.
    scoring_array = [[0 for j in range(len(s2) + 1)] for i in range(len(s1) + 1)]

    # fill in the first row and column with starting values to reflect indels (0, -1, -2, -3, ...).
    for i in range(1, len(s1) + 1):
        scoring_array[i][0] = i * indel
    for j in range(1, len(s2) + 1):
        scoring_array[0][j] = j * indel

    # fill in the rest of the matrix.
    for i in range(1, len(s1) + 1):
        for j in range(1, len(s2) + 1):
            
            # calculate the possible scores reflecting each possible alignment (a match, mismatch, insertion, or deletion).
            match_score = scoring_array[i-1][j-1] + S(s1[i-1], s2[j-1])
            insertion_score = scoring_array[i][j-1] + indel
            deletion_score = scoring_array[i-1][j] + indel
            
            val = max(match_score, deletion_score, insertion_score)

            scoring_array[i][j] = val
    
    # initialize the final alignment strings.
    alignment1, alignment2 = '', ''
    i, j = len(s1), len(s2)
    while (i > 0 or j > 0):
        # check for the score of the four possible alignment choices.

        # match or mismatch between the sequences.
        if (i > 0 and j > 0 and scoring_array[i][j] == scoring_array[i-1][j-1] + S(s1[i-1], s2[j-1])):
            # align the characters in both sequences.
            alignment1 += s1[i-1]
            alignment2 += s2[j-1]
            i -= 1
            j -= 1

        # insertion in the first sequence.
        elif (i > 0 and scoring_array[i][j] == scoring_array[i-1][j] + indel):
            # align a gap in the second sequence.
            alignment1 += s1[i-1]
            alignment2 += '-'
            i -= 1

        # deletion in the first sequence.
        elif (j > 0 and scoring_array[i][j] == scoring_array[i][j-1] + indel):
            # align a gap in the first sequence.
            alignment1 += '-'
            alignment2 += s2[j-1]
            j -= 1

    # reverse the strings, as we added from the end of both sequences.
    alignment1 = alignment1[::-1]
    alignment2 = alignment2[::-1]

    # return the optimal alignment and its score.
    return alignment1, alignment2, scoring_array[len(s1)][len(s2)]

if __name__ == '__main__':
    # get two sequences from command line arguments.
    s1 = sys.argv[1]
    s2 = sys.argv[2]

    # align and print.
    print(align(s1, s2))

    # thank you professor! see you in 3500 next semester




