"""
Generator for Rocket Nozzle Contours
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math
from scipy.optimize import fsolve
import time


def conical_nozzle(r_t, r_ch, eta, alpha_div=15, alpha_con=25, rf_c=1, rf_d=1, n_steps=1000):
    """Generates the contour of a conical nozzle.

    :param r_t: radius of throat
    :param r_ch: radius of chamber
    :param eta: area expansion factor
    :param alpha_div: divergent half-angle (degrees)
    :param alpha_con: convergent half-angle (degrees)
    :param rf_c: radius factor convergent side
    :param rf_d: radius factor divergent side
    :param n_steps: number of steps
    :return: Pandas dataFrame of x-y coordinates
    """
    alpha_div = math.radians(alpha_div)
    alpha_con = math.radians(alpha_con)
    r_e = eta**0.5 * r_t
    nozzle = pd.DataFrame(np.zeros((n_steps, 3)), columns={'x', 'y', 'section'})
    nozzle.loc[0, 'y'] = r_ch
    l_1 = (r_ch - r_t * (1 + rf_c * (1 - math.cos(alpha_con)))) / math.tan(alpha_con)
    l_2 = r_t * rf_c * math.sin(alpha_con)
    l_3 = r_t * rf_d * math.sin(alpha_div)
    l_4 = (r_e - r_t * (1 + rf_d * (1 - math.cos(alpha_div)))) / math.tan(alpha_div)
    l_sum = l_1 + l_2 + l_3 + l_4
    print(((r_t * (math.sqrt(eta) - 1) + r_t * rf_d * (1/math.cos(math.radians(15)) - 1)) / (math.tan(math.radians(15)))))
    print(l_3+l_4)

    # for index, row in nozzle.iterrows():
    #     row['x'] = row.name * (l_sum / n_steps)
    #
    # def cont_iter(x):
    #     if x <= l_1:
    #         return r_ch - nozzle['x'] * math.tan(alpha_con), 0
    #     if l_1 < x <= l_1 + l_2:
    #         return r_t * (1 + rf_c * (1 - math.cos(math.asin((l_2 - nozzle['x'] + l_1) / (r_t * rf_c))))), 1
    #     if l_1 + l_2 < x <= l_1 + l_2 + l_3:
    #         return r_t * (1 + rf_d * (1 - math.cos(math.asin((nozzle['x'] - l_1 - l_2) / (r_t * rf_d))))), 2
    #     if l_1 + l_2 + l_3 < x:
    #         return r_t * (1 + rf_d * (1 - math.cos(alpha_div))) + (nozzle['x'] - (l_1 + l_2 + l_3)) * math.tan(
    #             alpha_div), 3
    # nozzle['y'] = cont_iter(nozzle['x'])

    for index, row in nozzle.iterrows():
        row['x'] = row.name * (l_sum / n_steps)
        if row['x'] <= l_1:
            row['section'] = 0
            row['y'] = r_ch - row['x'] * math.tan(alpha_con)
        if l_1 < row['x'] <= l_1 + l_2:
            row['section'] = 1
            row['y'] = r_t * (1 + rf_c * (1 - math.cos(math.asin((l_2 - row['x'] + l_1) / (r_t * rf_c)))))
        if l_1 + l_2 < row['x'] <= l_1 + l_2 + l_3:
            row['section'] = 2
            row['y'] = r_t * (1 + rf_d * (1 - math.cos(math.asin((row['x'] - l_1 - l_2) / (r_t * rf_d)))))
        if l_1 + l_2 + l_3 < row['x']:
            row['section'] = 3
            row['y'] = r_t * (1 + rf_d * (1 - math.cos(alpha_div))) + (row['x'] - (l_1 + l_2 + l_3)) * math.tan(alpha_div)

    return nozzle


def bell_nozzle(r_t, r_ch, eta, alpha_div=15, alpha_con=25, alpha_end=12, rf_c=1.5, rf_d=0.382, len_per=0.8, n_steps=1000):
    alpha_div = math.radians(alpha_div)
    alpha_con = math.radians(alpha_con)
    alpha_end = math.radians(alpha_end)
    r_e = eta ** 0.5 * r_t
    nozzle = pd.DataFrame(np.zeros((n_steps, 3)), columns={'x', 'y', 'section'})
    nozzle.loc[0, 'y'] = r_ch
    l_1 = (r_ch - r_t * (1 + rf_c * (1 - math.cos(alpha_con)))) / math.tan(alpha_con)
    l_2 = r_t * rf_c * math.sin(alpha_con)
    l_3 = r_t * rf_d * math.sin(alpha_div)
    l_4 = (((r_t * (math.sqrt(eta) - 1) + r_t * rf_d * (1/math.cos(math.radians(15)) - 1)) / (math.tan(math.radians(15))
                                                                                              )) - l_3) * len_per
    l_sum = l_1 + l_2 + l_3 + l_4
    print(l_1, l_1 + l_2, l_1 + l_2 + l_3, l_1 + l_2 + l_3 + l_4)

    def equations(p):
        a, b, c = p
        return (b * math.sqrt((l_1 + l_2 + l_3) - a) + c - (r_t * (1 + rf_d * (1 - math.cos(alpha_div)))),
                b / (2 * math.sqrt((l_1 + l_2 + l_3) - a)) - math.tan(alpha_div),
                b / (2 * math.sqrt(l_sum - a)) - math.tan(alpha_end))

    a, b, c = fsolve(equations, (-10, -10, -10))
    print(a, b, c)

    for index, row in nozzle.iterrows():
        row['x'] = row.name * (l_sum / n_steps)
        if row['x'] <= l_1:
            row['section'] = 0
            row['y'] = r_ch - row['x'] * math.tan(alpha_con)
        if l_1 < row['x'] <= l_1 + l_2:
            row['section'] = 1
            row['y'] = r_t * (1 + rf_c * (1 - math.cos(math.asin((l_2 - row['x'] + l_1) / (r_t * rf_c)))))
        if l_1 + l_2 < row['x'] <= l_1 + l_2 + l_3:
            row['section'] = 2
            row['y'] = r_t * (1 + rf_d * (1 - math.cos(math.asin((row['x'] - l_1 - l_2) / (r_t * rf_d)))))
        if l_1 + l_2 + l_3 < row['x']:
            row['section'] = 3
            row['y'] = b * math.sqrt(row['x'] - a) + c
    return nozzle


if __name__ == '__main__':
    print('This Module is not supposed to be used on it\'s own')
    start_time = time.time()
    example = conical_nozzle(0.01, 0.02, 2.5, n_steps=10000)
    print("time elapsed: {:.2f}s".format(time.time() - start_time))
    print(example)
    fig, ax = plt.subplots()
    ax.plot(example['x'], example['y'])
    ax.set(ylim=(0, 0.03))
    ax.set_aspect('equal')
    ax2 = ax.twinx()
    ax2.plot(example['x'], example['section'], 'r:', lw=0.5)
    ax.grid(True)
    plt.show()
