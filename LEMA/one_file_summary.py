# Summarizes all the steps of the simulation code in a single file

from   scipy.integrate import solve_ivp
import math
import matplotlib.pyplot as plt

# Coupling factor between projectile and coil pair
# as a function of turns in coil and projectile position
def dk_dx(N, x): 
    x0 = 1e-2
    N0 = 400
    a1, a2     =  0.21651925, 0.04803981
    b1, b2     = -0.34392363, 0.07699456
    k1, k2, k3 =  0.37254388, 0.04442989, -0.01831165
    a    = a1 + a2 * (N/N0)
    b    = b1 + b2 * (N/N0)
    kmax = k1 + k2 * (N/N0) + k3 * (N/N0)**2
    return kmax * math.exp(-(a*x/x0)**2 - (b*x/x0)**4) * \
            (- 2 * x * (a/x0)**2 - 4 * x**3 * (b/x0)**4)

# Inductance of coil pair as a funciton of number of turns
def Ls(N):
    N0 = 400
    A, B, C  = 0.88747166e-3, 4.05174052e-3, 2.13676405e-3
    return A * (N/N0)**1.5 + B * (N/N0)**2 + C * (N/N0)**2.5

# Self inducatance and current of permanent magnet equivalent loop
def Lp(): return 5.4101e-08 # [H]
def Ip(): return 8737.6     # [A]

# Resistance of coil pair with N turns
def R(N):
    # Constants
    A_t  = 0.6         # [mm2]   Effective wire crossection area
    r_i  = 7.5         # [mm]    Inner radius of coil
    rho  = 61.02e-3    # [Ohm/m] Resistance of wire

    # Dependent variables
    d    = math.sqrt(A_t * N)               # [mm]  Thickness of coil
    r_o  = r_i+ d                         # [mm]  Outer radius of coil
    vol  = d * math.pi * (r_o**2 - r_i**2)  # [mm3] Volume of coil
    l    = 1e-3 * vol / A_t               # [m]   Length of wire

    return 2 * rho * l

def dy_dt(t, y, param):
    x, v, u, i = y
    m, C, N = param

    # ODE system describing accelerator
    dx_dt =   v
    dv_dt =   dk_dx(N, x) * math.sqrt(Lp()*Ls(N)) * Ip() * i / m
    du_dt = - i / C
    di_dt = - Ip() * dk_dx(N, x) * math.sqrt(Lp()*Ls(N)) * v / Ls(N) \
            + u / Ls(N) \
            - R(N) * i / Ls(N)

    return (dx_dt, dv_dt, du_dt, di_dt)

# Contstants
m  = 118.2e-3   # [kg] Projectile mass
C  = 560e-6     # [F]  Capacitor capacitance

# Parameters
N  =   700      #       Turns of wire in coils
dx =  -4e-2     # [m]   Activation position
u0 =   300      # [V]   Starting capacitor voltage
v0 =   0        # [m/s] Starting projectile speed

# Solve ODE system
y0  = (dx, v0, u0, 0)  # Initial values
tspan = (0, 0.05)      # Timespan
param = m, C, N
result = solve_ivp(lambda t, y: dy_dt(t, y, param), tspan, y0, rtol=1e-8)

# Plot result
t = result.t
(x, v, u, i) = result.y
plt.subplot(2,1,1)
plt.plot(t, v)
plt.subplot(2,1,2)
plt.plot(t, i)
plt.show()