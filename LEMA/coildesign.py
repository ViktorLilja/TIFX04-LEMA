import math

n = 670     # Number of turns
A_c = 0.6   # [mm2] Effective wire area
r   = 7.5   # [mm] Inner radius of coil


d = math.sqrt(A_c * n)  # [mm] Thickness of coil
R = r + d               # [mm] Outer radius of coil
l = 1e-3 * d * math.pi * (R**2 - r**2) / A_c    # [m] Length of wire

print("Number of turns %d" % n)
print("Coil thickness %.3f mm" % d)
print("Length of wire %.3f m" % l)