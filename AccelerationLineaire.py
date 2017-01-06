import numpy as np
from scipy.integrate import quad
from scipy.optimize import ridder

g = 9.81


def log(s):
    print('    {0}'.format(s))


def integ_table(table, t, res_0):
    res = [res_0]
    for i in range(1, len(t)):
        dt = t[i] - t[i - 1]
        res.append(res[i - 1] + table[i] * dt)
    return res


def etape_1(theta_final, a_retour, acc_ang_max):
    delta_theta = np.arctan(a_retour / g)
    T = 2 * np.sqrt((theta_final + delta_theta) / acc_ang_max)
    T_demi = T / 2

    f_1 = lambda x: g * np.tan(theta_final - 0.5 * acc_ang_max * x ** 2)
    f_2 = lambda x: g * np.tan(0.5 * acc_ang_max * (T - x) ** 2 - delta_theta)

    ts = np.linspace(0, T, 200)

    def a(t):
        if t > T_demi:
            return f_2(t)
        return f_1(t)

    les_a = [a(t) for t in ts]
    les_v = integ_table(les_a, ts, 0)
    les_x = integ_table(les_v, ts, 0)

    def theta(t):
        if t < T_demi:
            return 0.5 * acc_ang_max * t ** 2
        return - 0.5 * acc_ang_max * (T - t) ** 2 + theta_final + delta_theta

    thetas = [theta(t) for t in np.linspace(0, T, 200)]

    res = [[les_x[i], les_v[i], les_a[i]] for i in range(len(les_a))]

    return res, ts, thetas


def etape_2(a_r, theta, v_1, x_1=0, t0=0):
    T = v_1 / a_r
    ts = np.linspace(t0, t0 + T, 400)
    thetas = [theta for _ in ts]
    return [(x_1 + v_1 * t - 0.5 * a_r * t ** 2, v_1 - a_r * t, -a_r) for t in np.linspace(0, T, 400)], ts, thetas


def etape_3(a_retour, theta_final, acc_ang_max, v_e2, x_e2, t0):
    delta_theta = np.arctan(a_retour / g)
    theta_1, theta_2 = theta_final + delta_theta, theta_final
    delta = theta_2 - theta_1
    T = 2 * np.sqrt(abs(delta) / acc_ang_max)
    T_demi = T / 2
    sg = delta / abs(delta)

    ts = np.linspace(0, T, 200)

    f_1 = lambda x: g * np.tan(theta_final - (sg * 0.5 * acc_ang_max * (x ** 2) + theta_1))
    f_2 = lambda x: g * np.tan(theta_final - (-sg * 0.5 * acc_ang_max * (T - x) ** 2 + theta_2))

    def a(t):
        if t > T_demi:
            return f_2(t)
        return f_1(t)

    def theta(t):
        if t < T_demi:
            return sg * 0.5 * acc_ang_max * t ** 2 + theta_1
        return -sg * 0.5 * acc_ang_max * (T - t) ** 2 + theta_2

    thetas = [theta(t) for t in np.linspace(0, T, 200)]

    les_a = [a(t) for t in ts]
    les_v = integ_table(les_a, ts, v_e2)
    les_x = integ_table(les_v, ts, x_e2)
    res = [[les_x[i], les_v[i], les_a[i]] for i in range(len(les_a))]

    t = [e + t0 for e in np.linspace(0, T, 200)]
    return res, t, thetas


def etape_5(a_r, theta_final, acc_ang_max, delta_v):

    delta_theta = np.arctan(a_r / g)

    def vT(theta, n=100):
        theta_1, theta_2 = theta_final, theta_final - theta
        delta = theta_2 - theta_1
        T = 2 * np.sqrt(abs(delta) / acc_ang_max)
        T_demi = T / 2
        sg = 1
        if delta<0:
            sg = -1

        f_1 = lambda x: g * np.tan(theta_final - (sg * 0.5 * acc_ang_max * (x ** 2) + theta_1))
        f_2 = lambda x: g * np.tan(theta_final - (-sg * 0.5 * acc_ang_max * (T - x) ** 2 + theta_2))

        v1, v2 = quad(f_1, 0, T_demi, limit=n // 2), quad(f_2, T_demi, T, limit=n // 2)
        return v1[0] + v2[0]

    def make_x_v_a(theta1, theta2, v0):

        theta_1, theta_2 = theta_final + theta1, theta_final + theta2
        delta = theta_2 - theta_1
        T = 2 * np.sqrt(abs(delta) / acc_ang_max)
        T_demi = T / 2
        sg = delta / abs(delta)

        f_1 = lambda x: g * np.tan(theta_final - (sg * 0.5 * acc_ang_max * (x ** 2) + theta_1))
        f_2 = lambda x: g * np.tan(theta_final - (-sg * 0.5 * acc_ang_max * (T - x) ** 2 + theta_2))

        g_1 = lambda t: quad(f_1, 0, t)[0] + v0
        g_2 = lambda t: quad(f_2, T_demi, t)[0]

        x_1 = lambda u: quad(g_1, 0, u)[0]
        x_2 = lambda u: quad(g_2, T_demi, u)[0] + (u - T_demi) * g_1(T_demi)

        def x(t):
            if t > T_demi:
                return x_1(T_demi) + x_2(t)
            return x_1(t)

        def v(t):
            if t > T_demi:
                return g_1(T_demi) + g_2(t)
            return g_1(t)

        def a(t):
            if t > T_demi:
                return f_2(t)
            return f_1(t)

        def thetas(t):
            if t < T_demi:
                return sg * 0.5 * acc_ang_max * t ** 2 + theta_1
            return -sg * 0.5 * acc_ang_max * (T - t) ** 2 + theta_2

        return x, v, a, thetas, T, T_demi

    log('vT(0) = {0}'.format(vT(0)))
    log('vT(delta_theta) = {0}'.format(vT(delta_theta)))
    log('delta_v = {0}'.format(delta_v))

    if vT(delta_theta) < -delta_v / 2:
        log('Will go easy way')
        v_T = vT(delta_theta)
        x1, v1, a1, thetas1, T_1, T_demi_1 = make_x_v_a(0, -delta_theta, delta_v)
        x2, v2, a2, thetas2, T_2, T_demi_2 = make_x_v_a(-delta_theta, 0, -v_T)



        return None

    log('Will go hard way')
    log('Looking for theta value')

    f = lambda x: vT(x) + (delta_v / 2)

    log('f(0) = {0}'.format(f(0)))
    log('f(delta_theta) = {0}'.format(f(delta_theta)))
    log('delta_v/2 = {0}'.format(delta_v / 2))

    sol = ridder(f, 0, delta_theta)

    log('sol = {0}'.format(sol))
    log('| vT(sol) - |delta_theta| | = {0}'.format(abs(vT(sol, 1000) + delta_v/2)))

    x1, v1, a1, thetas1, T_1, T_demi_1 = make_x_v_a(0, -sol, delta_v)
    x2, v2, a2, thetas2, T_2, T_demi_2 = make_x_v_a(-sol, 0, delta_v / 2)

    x1_max = x1(T_1)

    t_1 = np.linspace(0, T_1, 100)
    res_1, thetas_1 = [(x1(t) - (x2(T_2) + x1_max), v1(t), a1(t)) for t in t_1], [thetas1(t) for t in t_1]

    t_2 = np.linspace(0, T_2, 100)
    res_2, thetas_2 = [(x2(t) - x2(T_2), v2(t), a2(t)) for t in t_2], [thetas2(t) for t in t_1]

    res = res_1 + res_2
    t_1, t_2 = [e for e in t_1], [e + T_1 for e in t_2]
    t = t_1 + t_2
    thetas = thetas_1 + thetas_2

    return res, t, thetas


def etape_4(x3, x5, v3, t0, theta_f):
    T = (x5 - x3) / v3
    t = np.linspace(0, T, 200)
    return [(x3 + v3 * t, v3, 0) for t in t], t + t0, [theta_f for _ in t]




