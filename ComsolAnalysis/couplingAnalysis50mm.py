import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import minimize

# --- READ DATA --- #

data = pd.read_csv("ComsolAnalysis/data/20220325_Coupling50mm.csv")

# Coordinates (projectile pos, coil turns)
x = np.array(data["pos"])
n = np.array(list(range(100, 801, 100)))
n_m, x_m = np.meshgrid(n, x)

# Total magnetic energy for configuration
n_str = list(map(lambda x : str(x), n))
W = np.array(data[n_str])


# --- CALCULATE COUPLING COEFFICIENT FROM DATA --- #

# Magnet energy (from Comsol)
Wm = 2.0652 # J

# Magnet surface current (from Comsol)
Im = 8737.6 # A
print("Magnet Current:", Im, "A")

# Magnet equivalent inductance
Lm = 2 * Wm / Im**2
print("Magnet inductance:", Lm, "H")
print(Lm)

# Coil current (from Comsol)
Ic = 20 * (500/n_m)

# Coil inductance (from coilAnalysis)
n0 = 400
A  = 0.87887e-3
B  = 4.01249e-3
C  = 2.11607e-3
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
    (a1, a2, b1, b2, k0, k1, k2, k3) = param
    a    = a1 + a2 * (n/n0)
    b    = b1 + b2 * (n/n0)
    kmax = k1 + k2 * (n/n0) + k3 * (n/n0)**2
    return k0 + kmax * np.exp(- (a*x/x0)**2 - (b*x/x0)**4)
        

# Optimize
def error(param):
    return np.sum((k - k_fit(x_m, n_m, param))**2)

res = minimize(error, (1, 1, 1, 1, 1, 1, 1, 1), tol=1e-9, options={"maxiter": 10000})
param = res.x

# Plot all coupling curves
plt.plot(1e2*x[3:], k[3:,1], "rx")
plt.plot(1e2*x[3:], k[3:,3], "gx")
plt.plot(1e2*x[3:], k[3:,5], "bx")
plt.plot(1e2*x[3:], k_fit(x_m, n_m, param)[3:,1], "r--")
plt.plot(1e2*x[3:], k_fit(x_m, n_m, param)[3:,3], "g--")
plt.plot(1e2*x[3:], k_fit(x_m, n_m, param)[3:,5], "b--")
plt.ylabel("Kopplingskoefficient, k")
plt.xlabel("Projektilposition, x (cm)")
plt.legend(
    list(map(lambda x: "Beräknat i COMSOL, " + x + " varv ", n_str[1:6:2])) + 
    list(map(lambda x: "Kurvanpassning, " + x + " varv ", n_str[1:6:2]))
    )
plt.grid()
plt.show()

# Write to csv for plotting
with open("ComsolAnalysis/output/COMSOL_k.csv", "w") as f:
    # Header
    f.write("x (cm),data 200,data 400,data 600,fit 200, fit 400, fit 600\n")


    for i in range(3, len(x)):
        f.write("%e," % x[i])
        f.write("%e," % k[i,1])
        f.write("%e," % k[i,3])
        f.write("%e," % k[i,5])
        f.write("%e," % k_fit(x_m, n_m, param)[i,1])
        f.write("%e," % k_fit(x_m, n_m, param)[i,3])
        f.write("%e\n" % k_fit(x_m, n_m, param)[i,5])

# kmax vs n
#plt.plot(n, k[-1,:], "o-")
#plt.plot(n, k_fit(0, n, param), "k--")
#plt.xlabel("Number of turns")
#plt.ylabel("Maximum coupling")
#plt.show()

# Sample curve
#i = 2
#plt.plot(x, k[:, i], "ro")
#plt.plot(x, k_fit(x, n[i], param), "k--")
#plt.title("%s turns" % n[i])
#plt.show()

#(a1, a2, b1, b2, k0, k1, k2, k3) = param
#print("x0 = %.8f" % x0 )
#print("n0 = %d" % n0 )
#print("a1 = %.8f" % a1 )
#print("a2 = %.8f" % a2 )
#print("b1 = %.8f" % b1 )
#print("b2 = %.8f" % b2 )
#print("k1 = %.8f" % k1 )
#print("k2 = %.8f" % k2 )
#print("k3 = %.8f" % k3 )
#print("a    = a1 + a2 * (n/n0)")
#print("b    = b1 + b2 * (n/n0)")
#print("kmax = k1 + k2 * (n/n0) + k3 * (n/n0)**2")
#print("k = kmax * np.exp(- (a*x/x0)**2 - (b*x/x0)**4)")
#print("dkdx = kmax * np.exp(-(a*x/x0)**2 - (b*x/x0)**4) * (-2*x*(a/x0)**2 -4*x**3*(b/x0)**4)")