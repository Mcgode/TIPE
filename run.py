import AccelerationLineaire as a_l
import numpy as np
import matplotlib.pyplot as plt
import Utils

from time import time

g = 9.81
al = 0.3 * g

max_acc_ang = 1
a_r = 1

theta_f = np.arccos(g / np.sqrt(al ** 2 + g ** 2))

T0 = time()

t_1 = time()
res1, t1, thetas_1 = a_l.etape_1(theta_f, a_r, max_acc_ang)

t_2 = time()
res2, t2, thetas_2 = a_l.etape_2(a_r, thetas_1[-1], res1[-1][1], res1[-1][0], t1[-1])

t_3 = time()
res3, t3, thetas_3 = a_l.etape_3(a_r, theta_f, max_acc_ang, res2[-1][1], res2[-1][0], t2[-1])

t_5 = time()
res5, t5, thetas_5, T = a_l.etape_5(a_r, theta_f, max_acc_ang, res3[-1][1], t3[-1])

if T is not None:
    t_s = np.linspace(0, T, 100)
    res_plus = [(res2[-1][0] - 0.5 * a_r * t ** 2, -a_r * t, -a_r) for t in t_s]
    delta_x, delta_v = -0.5 * a_r * T ** 2, -a_r * T
    for i in range(len(res3)):
        last_x, dt = res_plus[-1][0], 0
        if i > 0:
            dt = t3[i] - t3[i - 1]
            last_x = res3[i - 1][0]
        res3[i] = (last_x + (res3[i][1] + delta_v) * dt, res3[i][1] + delta_v, res3[i][2])
    res3 = res_plus + res3
    t3 = [t2[-1] + t for t in t_s] + [t + T for t in t3]
    thetas_3 = [thetas_3[0] for _ in t_s] + thetas_3

t_4 = time()
res4, t4, thetas_4 = a_l.etape_4(res3[-1][0], res5[0][0], res3[-1][1], t3[-1], theta_f)

T1 = time()

print('Calculation process took {0} seconds'.format(T1 - T0))
print('   Step 1 took {0} seconds'.format(t_2 - t_1))
print('   Step 2 took {0} seconds'.format(t_3 - t_2))
print('   Step 3 took {0} seconds'.format(t_5 - t_3))
print('   Step 4 took {0} seconds'.format(T1 - t_4))
print('   Step 5 took {0} seconds'.format(t_4 - t_5))

T2 = time()

t1, t2, t3, t4 = [e for e in t1], [e for e in t2], [e for e in t3], [e for e in t4]
t5 = [e + t4[-1] for e in t5]
t_s = t1 + t2 + t3 + t4 + t5

res = res1 + res2 + res3 + res4 + res5
x, v, a = [e[0] for e in res], [e[1] for e in res], [e[2] for e in res]

thetas = (180 / np.pi) * np.array(thetas_1 + thetas_2 + thetas_3 + thetas_4 + thetas_5)

T3 = time()

print('Data regroup process took {0} seconds'.format(T3 - T2))

T4 = time()
Utils.save_pos_and_theta(t_s, res, thetas, 500)
T5 = time()

print('Data save process took {0} seconds'.format(T5 - T4))

plt.plot(t_s, x, label='x(t)')
plt.plot(t_s, v, label='v(t)')
plt.plot(t_s, a, label='a(t)')

maxima = Utils.max_spe(res)
minima = Utils.min_spe(res)

plt.plot([t1[-1], t1[-1]], [np.floor(max(maxima)) + 1, np.ceil(min(minima)) - 1], color='cyan')
plt.plot([t2[-1], t2[-1]], [np.floor(max(maxima)) + 1, np.ceil(min(minima)) - 1], color='cyan')
plt.plot([t3[-1], t3[-1]], [np.floor(max(maxima)) + 1, np.ceil(min(minima)) - 1], color='cyan')
plt.plot([t4[-1], t4[-1]], [np.floor(max(maxima)) + 1, np.ceil(min(minima)) - 1], color='cyan')

plt.legend(loc=2)
plt.xlabel('Temps t en s')

plt.grid()
plt.show()

plt.plot(t_s, [(180 / np.pi) * theta_f for _ in t_s], label='theta_f', color='green')
plt.plot(t_s, thetas, label='theta(t)', color='red')

plt.plot([t1[-1], t1[-1]], [np.floor(np.max(thetas) / 5) * 5 + 5, 0], color='cyan')
plt.plot([t2[-1], t2[-1]], [np.floor(np.max(thetas) / 5) * 5 + 5, 0], color='cyan')
plt.plot([t3[-1], t3[-1]], [np.floor(np.max(thetas) / 5) * 5 + 5, 0], color='cyan')
plt.plot([t4[-1], t4[-1]], [np.floor(np.max(thetas) / 5) * 5 + 5, 0], color='cyan')

plt.legend(loc=5)
plt.xlabel('Temps t en s')
plt.grid()
plt.show()
