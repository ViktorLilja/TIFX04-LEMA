import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

# --- READ DATA --- #

data = pd.read_csv("ComsolAnalysis/data/20220325_Inductance.csv")

# Coordinates (projectile pos, coil turns)
gap = np.array(data["gap"])
n = np.array(list(range(100, 801, 100)))
n_m , gap_m = np.meshgrid(n, gap)

# Total magnetic energy for configuration
n_str = list(map(lambda x : str(x), n))
Wc = np.array(data[n_str])


# --- CALCULATE INDUCTANCE FROM DATA --- #

# Coil current (from Comsol parameter file)
Ic = 20 * (500/n_m)


# Coil inductance
Lc = 2 * Wc / Ic**2


# --- FIT CURVE TO INDUCTANCE --- #
n0 = 400
g0 = 1e-2
def L_fit(gap_m, n_m, param):
    (a, b, c, A, B, C) = param
    L0  = (A * (n_m/n0)**1.5 + B * (n_m/n0)**2 + C * (n_m/n0)**2.5) * 1e-3
    Cg   = 2 + a * (gap_m/g0-13)**1 + b * (gap_m/g0-13)**2 + c * (gap_m/g0-13)**3
    return Cg * L0


# Optimize
def error(param):
    return np.sum(((Lc - L_fit(gap_m, n_m, param))/Lc)**2)

res = minimize(error, (0, 1, 1, 1, 1, 1), tol=1e-9, options={"maxiter": 10000})
param = res.x


# Plot all data
plt.plot(1e3*gap, 1e3*Lc, "x-")
plt.plot(1e3*gap, 1e3*L_fit(gap_m, n_m, param), "k--")
plt.ylabel("Inductance (mH)")
plt.xlabel("Coil gap (mm)")
plt.legend(list(map(lambda x: "n = " + x, n_str)))
plt.show()

# k - curve
(a, b, c, A, B, C) = param
n_sample = 400
g_space = np.linspace(0.006, 0.016)
Cg   = 2 + a * (g_space/g0-13)**1 + b * (g_space/g0-13)**2 + c * (g_space/g0-13)**3
plt.plot(1e3*g_space, Cg)
plt.title("n = %d" % n_sample)
plt.ylabel("Gap factor")
plt.xlabel("Coil gap (mm)")
plt.show()


print("n0 = %d" % n0 )
print("g0 = %-8f" % g0 )
print("a  = %.8f" % a)
print("b  = %.8f" % b)
print("c  = %.8f" % c)
print("A  = %.8f" % A)
print("B  = %.8f" % B)
print("C  = %.8f" % C)
print("L0 = (A * (n/n0)**1.5 + B * (n/n0)**2 + C * (n/n0)**2.5) * 1e-3")
print("Cg = 2 + a * (gap/g0-13)**1 + b * (gap/g0-13)**2 + c * (gap/g0-13)**3")
print("L  = L0 * Cg")