import numpy as np

# Vacuum permeability
__mu0__ = 1.25663706212e-6

# Coil width and thickness [mm]
def width(n): return np.sqrt(0.6 * n)

# Length of wire in coil [m]
def length(n):
    d = width(n)
    return 5.26e-3 * d**2 * (d+15)

# Resistance of coil pair [Ohm]
def resistance(n): return 0.067 * 2 * length(n)

# Inductance of coil pair [H]
def inductance(gap, n): 
    n0 = 400
    g0 = 1e-2
    a  = -0.82546566
    b  = -0.14298914
    c  = -0.00683613
    A  = 0.28977120
    B  = 1.32294670
    C  = 0.69768164
    L0 = (A * (n/n0)**1.5 + B * (n/n0)**2 + C * (n/n0)**2.5) * 1e-3
    Cg = 2 + a * (gap/g0-13)**1 + b * (gap/g0-13)**2 + c * (gap/g0-13)**3
    return L0 * Cg

class Stage:

    # Initialize and calculate parameters
    def __init__(self, 
                 n,               # Number of turns per coil in stage
                 gap,             # Gap between coils 
                 x = 0,           # Position of coils along rail
                 C = 560e-6,      # Capacitance of capacitor
                 dx = -5e-2,      # Trigger position of stage
                 uC0 = 300,       # Initial capacitor voltage
                 i0 = 0,          # Initial coil current 
                 ):

        self.n      = n
        self.L      = inductance(gap, n)
        self.R      = resistance(n)
        self.x      = x
        self.C      = C
        self.dx     = dx
        self.uC0    = uC0
        self.i0     = i0

        self.hasTrigged = False

    # Returns true if coil if stage is activated given projectile position x
    def active(self, x):
        if self.hasTrigged: return True
        self.hasTrigged = x - self.x >= self.dx
        return self.hasTrigged
