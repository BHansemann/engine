'''
Generator for Rocket Nozzle Contours
'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math


def conical_nozzle(r_t, r_ch, eta, alpha_div=15, alpha_con=35, rf_t=1, n_steps=1000):
    alpha_div = math.radians(alpha_div)
    alpha_con = math.radians(alpha_con)
    r_e = (eta * r_t**2)**0.5
    nozzle = pd.DataFrame(np.zeros((n_steps, 3)), columns={'x', 'y', 'section'})
    nozzle.loc[0, 'y'] = r_ch
    l_1 = (r_ch - r_t * (1 + rf_t * (1 - math.cos(alpha_con)))) / math.tan(alpha_con)
    l_2 = r_t * rf_t * math.sin(alpha_con)
    l_3 = r_t * rf_t * math.sin(alpha_div)
    l_4 = (r_e - r_t * (1 + rf_t * (1 - math.cos(alpha_div)))) / math.tan(alpha_div)
    l_sum = l_1 + l_2 + l_3 + l_4
    for index, row in nozzle.iterrows():
        row['x'] = row.name * (l_sum / n_steps)
        row['section'] = 3
        row['x']
        if row['x'] <= l_1:
            row['section'] = 0
            row['y'] = r_ch - row['x'] * math.tan(alpha_con)
        if l_1 < row['x'] <= l_1 + l_2:
            row['section'] = 1
            row['y'] = r_t * (1 + (1 - math.cos(math.asin((l_2 - row['x'] + l_1) / (r_t * rf_t)))))
        if l_1 + l_2 < row['x'] <= l_1 + l_2 + l_3:
            row['section'] = 2
            row['y'] = r_t * (1 + (1 - math.cos(math.asin((row['x'] - l_1 - l_2) / (r_t * rf_t)))))
        if l_1 + l_2 + l_3 < row['x']:
            row['section'] = 3
            row['y'] = r_t * (1 + rf_t * (1 - math.cos(alpha_div))) + (row['x'] - (l_1 + l_2 + l_3)) * math.tan(alpha_div)

    return nozzle


if __name__ == '__main__':
    print('This Module is not supposed to be used on it\'s own')
    example = conical_nozzle(0.01, 0.02, 2.5)
    print(example)
    fig, ax = plt.subplots()
    ax.plot(example['x'], example['y'])
    ax.set(ylim=(0, 0.03))
    ax.set_aspect('equal')
    ax2 = ax.twinx()
    ax2.plot(example['x'], example['section'], 'r:', lw=0.5)
    ax.grid(True)
    plt.show()
