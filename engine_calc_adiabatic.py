# adiabatic approximation of thermodynamic properties inside the engine
import CoolProp.CoolProp as CP
import numpy as np
import pandas as pd

# inputs
p_ch = 100000       # chamber pressure in Pa
p_e = 10000         # exit pressure in Pa
n_ps = 100          # number of pressure steps
T_ch = 2500         # chamber temperature in Kelvin
F_thrust = 500      # Thrust in N
x_n2 = 0.34269      # Mole fraction of exhaust products, via NASA CEA
x_co = 0.23055
x_h2 = 0.21104
x_h2o = 0.17950
x_co2 = 0.03140
x_h = 0.00482

components_string = "n2[{}]&co[{}]&h2[{}]&h2o[{}]&co2[{}]&h[{}]".format(x_n2, x_co, x_h2, x_h2o, x_co2, x_h)
print(components_string)

p_steps = list(range(p_ch, p_e+1, (p_ch - p_e)//n_ps))

data = np.zeros((n_ps, 7))
print(data)