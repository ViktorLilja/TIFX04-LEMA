# Calculate coil parameters given numer of turns

import math

#    SP1 SP2 SP3 SP4 SP5
# n  668 470 358 313 274
# v0 0.0 3.1 5.6 7.0 8.5
n   = 668     # Number of turns
A_c = 0.6   # [mm2] Effective wire area
r   = 7.5   # [mm] Inner radius of coil


d = math.sqrt(A_c * n)  # [mm] Thickness of coil
R = r + d               # [mm] Outer radius of coil
l = 1e-3 * d * math.pi * (R**2 - r**2) / A_c    # [m] Length of wire

resistance = 0.06102 * l

print("Number of turns %d" % n)
print("Coil thickness %.3f mm" % d)
print("Length of wire %.3f m" % l)
print("Estimated resistance: %.3f Ohm" % resistance)