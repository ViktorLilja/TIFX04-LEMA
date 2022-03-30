import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

# --- READ DATA --- #

data = pd.read_csv("ComsolAnalysis/data/20220325_Coupling2x20mm.csv")

# Coordinates (projectile pos, coil turns)
x = np.array(data["pos"])
n = np.array(list(range(100, 801, 100)))
n_m, x_m = np.meshgrid(n, x)

# Total magnetic energy for configuration
n_str = list(map(lambda x : str(x), n))
W = np.array(data[n_str])


# --- CALCULATE COUPLING COEFFICIENT FROM DATA --- #

# Magnet energy (from Comsol)
Wm = 2.9387 # J

# Magnet surface current (from Comsol)
Im = 10345 # A
print("Magnet Current:", Im, "A")

# Magnet equivalent inductance
Lm = 2 * Wm / Im**2
print("Magnet inductance:", Lm, "H")

# Coil current (from Comsol)
Ic = 20 * (500/n_m)

# Coil inductance (from coilAnalysis)
n0 = 400
A  = 0.87084e-3
B  = 3.97581e-3
C  = 2.09672e-3
Lc  = A*(n_m/n0)**1.5 + B*(n_m/n0)**2 + C*(n_m/n0)**2.5

# Magnetic energy form coil
Wc = 0.5 * Lc * Ic**2

# Mutual inductance
M = (W - Wm - Wc) / (Im * Ic)

# Coupling factor
k = M / np.sqrt(Lc * Lm)

# --- FIT CURVE TO COUPLING COEFFICIENT --- #
x0 = 1e-2
n0 = 400
def k_fit(x, n, param):
    (a1, a2, a3,  b1, b2, b3, c1, c2, c3) = param
    a = a1 + a2 * (n/n0) + a3 * (n/n0)**2
    b = b1 + b2 * (n/n0) + b3 * (n/n0)**2
    c = c1 + c2 * (n/n0) + c3 * (n/n0)**2
    return a * np.sin(c*x/x0/b) * np.exp(- (x/x0/b)**2)

# Optimize
def error(param):
    return np.sum((k - k_fit(x_m, n_m, param))**2)

res = minimize(error, (1, 1, 1, 1, 1, 1, 1, 1, 1), tol=1e-9, options={"maxiter": 10000})
param = res.x

# Plot all coupling curves
plt.plot(1e2*x, k, "o-")
plt.plot(1e2*x, k_fit(x_m, n_m, param), "k--")
plt.ylabel("Coupling coefficent")
plt.xlabel("Projectile position (cm)")
plt.legend(list(map(lambda x: "n = " + x, n_str)))
plt.show()

# Sample curve
i = 4
plt.plot(x, k[:, i], "ro")
plt.plot(x, k_fit(x, n[i], param), "k--")
plt.title("%s turns" % n[i])
plt.show()

(a1, a2, a3,  b1, b2, b3, c1, c2, c3) = param
print("x0 = %.8f" % x0 )
print("n0 = %d" % n0 )
print("a1 = %.8f" % a1 )
print("a2 = %.8f" % a2 )
print("a3 = %.8f" % a3 )
print("b1 = %.8f" % b1 )
print("b2 = %.8f" % b2 )
print("b3 = %.8f" % b3 )
print("c1 = %.8f" % c1 )
print("c2 = %.8f" % c2 )
print("c3 = %.8f" % c3 )
print("a = a1 + a2 * (n/n0) + a3 * (n/n0)**2")
print("b = b1 + b2 * (n/n0) + b3 * (n/n0)**2")
print("c = c1 + c2 * (n/n0) + c3 * (n/n0)**2")
print("k = a * np.sin(c*x/x0/b) * np.exp(-(x/x0/b)**2)")
print("dkdx = (a*c/x0/b) * (np.cos(c*x/x0/b) - (2*x/x0/b) * np.sin(c*x/x0/b)) * np.exp(-(x/x0/b)**2)")