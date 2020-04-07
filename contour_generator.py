'''
Generator for Rocket Nozzle Contours
'''
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import math


def conical_nozzle(r_t, r_ch, eta, alpha_div=15, alpha_con=20, rf_t=1, n_steps=100):
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
        if row['x'] < l_1 + l_2 + l_3:
            row['section'] = 2
        if row['x'] < l_1 + l_2:
            row['section'] = 1
        if row['x'] < l_1:
            row['section'] = 0

    for index, row in nozzle.iterrows():
        if row['section'] == 0:
            row['y'] = r_ch - row['x'] * math.sin(alpha_con)

    return nozzle


if __name__ == '__main__':
    print('This Module is not supposed to be used on it\'s own')
    example = conical_nozzle(0.01, 0.02, 5)
    print(example)
    plt.plot(example['x'], example['y'])
    plt.plot(example['x'], example['section'])
    plt.show()
