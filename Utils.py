
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


def save_all_data(t, res, thetas, filename='/Users/max/Desktop/data.csv'):

    f = open(filename, 'w')
    data = [[t[i], res[i][0], res[i][1], res[i][2], thetas[i]] for i in range(len(t))]
    for e in data:
        string = ''
        for p in e:
            string += ("\"%.14f\"," % p).replace('.', ',');
        string += '\n'
        f.write(string)
    f.close()
