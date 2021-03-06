from numpy import pi


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


def save_all_data(t, res, thetas, smoothing=0, filename='/Users/max/Desktop/data.csv'):

    f = open(filename, 'w')
    f.write('t,x,v,a,theta,\n')
    data = [[t[i], res[i][0], res[i][1], res[i][2], thetas[i]] for i in range(len(t))]
    if smoothing > 0:
        data = smooth(data, smoothing)
    for e in data:
        string = ''
        for p in e:
            string += ("\"%.14f\"," % p).replace('.', ',')
        string += '\n'
        f.write(string)
    f.close()


def smooth(data, n):
    assert n > 0
    step = (data[-1][0] - data[0][0]) / n
    i = 0
    result = [data[0]]
    for j in range(1, n - 1):
        while step * j > data[i][0]:
            i += 1
        d1, d2 = data[i - 1], data[i]
        dt = d2[0] - d1[0]
        f = (j * step - d1[0]) / dt
        r = [d2[k] + f * (d2[k] - d1[k]) for k in range(len(d1))]
        result.append(r)
    result.append(data[-1])
    return result


def save_pos_and_theta(t, res, thetas, n=0):

    file = open('/Users/max/Desktop/res.csv', 'w')
    print(len(res), len(t), len(thetas))
    data = [[t[i], res[i][0], (thetas[i] * pi / 180)] for i in range(len(t))]
    maximum = data[0][1]
    for e in data:
        if e[1] > maximum:
            maximum = e[1]
    print("Maximum: %.2f" % maximum)
    if n > 0:
        data = smooth(data, n)
    for e in data:
        string = ("\"%.14f\",\"%.14f\"\n" % (e[0], e[1])).replace('.', ',')
        file.write(string)
    file.close()


def save_speed(res, t, n=200):
    data = smooth([[t[i], res[i][1]] for i in range(len(res))], n)
    with open('/Users/max/Desktop/speed.txt', 'w') as f:
        f.write("{0}\n".format(int(1000 * data[1][0])))
        for e in data:
            f.write('{0}\n'.format(int(1000 * e[1])))
        f.close()



