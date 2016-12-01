

def max_spe(t):
    l = len(t[0])
    maxima = [e for e in t[0]]

    for i in range(len(t)):
        for j in range(l):
            if t[i][j] > maxima[j]:
                maxima[j] = t[i][j]

    return maxima

def min_spe(t):
    l = len(t[0])
    minima = [e for e in t[0]]

    for i in range(len(t)):
        for j in range(l):
            if t[i][j] < minima[j]:
                minima[j] = t[i][j]

    return minima
