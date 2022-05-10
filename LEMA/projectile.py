# Projectile class, contains parameters relating to a specific
# type of permanent magnet projectile

import numpy as np

# Vacuum permeability
__mu0__ = 1.25663706212e-6

class Projectile:

    def __init__(self,
                 type,                  # 50mm or 2x20mm
                 x0    = -1e-2,         # Projectile starting position
                 xdot0 = 0,             # Projectile starting speed
                 ):
 
        self.type   = type
        self.x0     = x0
        self.xdot0  = xdot0

        if self.type=="50mm":
            self.m  = 51.6e-3 + 68.4e-3   # Projectile mass
            self.L  = 5.4101e-08          # Equivalent coil inductance
            self.I  = 8737.6              # Equivalent coil current
        elif self.type=="2x20mm":
            self.m  = 51.6e-3 + 60.8e-3
            self.L  = 5.4919e-08
            self.I  = 10345
        else:
            raise Exception("Unknown projectile type.")
    
    # Get appropriate coil gap for this type of magnet
    def gap(self):
        if self.type=="50mm":
            return 13e-3
        elif self.type=="2x20mm":
            return 14e-3
        else:
            raise Exception("Unknown projectile type.")

    def dk_dx(self, n, x):
        if self.type=="50mm":
            x0 = 1e-2
            n0 = 400
            a1 = 0.21651925
            a2 = 0.04803981
            b1 = -0.34392363
            b2 = 0.07699456
            k1 = 0.37254388
            k2 = 0.04442989
            k3 = -0.01831165
            a    = a1 + a2 * (n/n0)
            b    = b1 + b2 * (n/n0)
            kmax = k1 + k2 * (n/n0) + k3 * (n/n0)**2
            return kmax * np.exp(-(a*x/x0)**2 - (b*x/x0)**4) * \
                   (- 2 * x * (a/x0)**2 - 4 * x**3 * (b/x0)**4)

        elif self.type=="2x20mm":
            x0 = 1e-2
            n0 = 400
            a1 = 0.34084244
            a2 = -0.19481632
            a3 = 0.04373308
            b1 = 2.19392915
            b2 = 1.47545206
            b3 = -0.26377058
            c1 = 2.37595129
            c2 = 0.85278299
            c3 = -0.23910952
            a = a1 + a2 * (n/n0) + a3 * (n/n0)**2
            b = b1 + b2 * (n/n0) + b3 * (n/n0)**2
            c = c1 + c2 * (n/n0) + c3 * (n/n0)**2
            return (a*c/x0/b) * \
                   (np.cos(c*x/x0/b) - (2*x/x0/b) * np.sin(c*x/x0/b)) * \
                   np.exp(-(x/x0/b)**2)
        
        else:
            raise Exception("Unknown projectile type.")