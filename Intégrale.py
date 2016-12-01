import scipy.integrate as intg
import numpy as np
import matplotlib.pyplot as plt


def f1(t):
    omega = 0.5
    return 1 / np.cos(0.5 * omega * t*t)


def f2(t):
    omega = 0.5
    return np.tan(0.5 * omega * t*t)

omega = 0.5
g = 9.81
result = []
n = 1000

for e in np.linspace(0, 1, n):
    alpha = e # alpha = aR / g
    T = (np.arcsin(alpha) * 2 / omega) ** 0.5
    r1, r2 = intg.quad(f1, 0, T)[0] * g * alpha, intg.quad(f2, 0, T)[0] * g
    rapport = r1 / r2
    result.append(rapport)


plt.plot(np.linspace(0, 1, n), result)
plt.show()
