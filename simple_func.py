def hamming(x1, x2, m):
    j = 0
    for i in range(len(x1)):
        if x1[i] != x2[i]:
            j = j + 1
            if j > m:
                return False
    return True

# return Hamming distance between two strings
def hamming_dist(x1, x2):
    j = 0
    for i in range(len(x1)):
        if x1[i] != x2[i]:
            j = j + 1
    return j
