# Prandtl-Meyer Expansion Tool
# TI-84 CE Python compatible

import math

def nu(M, gamma):
    a = math.sqrt((gamma+1)/(gamma-1))
    b = math.sqrt((gamma-1)/(gamma+1)*(M**2 - 1))
    return a*math.atan(b) - math.atan(math.sqrt(M**2 - 1))

def solve_M_from_nu(nu_target, gamma):
    M_low = 1.0001
    M_high = 20.0

    for i in range(60):
        M_mid = 0.5*(M_low + M_high)
        if nu(M_mid, gamma) > nu_target:
            M_high = M_mid
        else:
            M_low = M_mid

    return 0.5*(M_low + M_high)

# -------------------------
# USER INPUT
# -------------------------

gamma = float(input("Gamma (1.4 for air) = "))

print("Mode:")
print("1 = Input Mach")
print("2 = Input nu (deg)")
mode = int(input("Mode #: "))

if mode == 1:
    M = float(input("Mach number M = "))
    nu_val = nu(M, gamma)
    print("nu (deg) = ", nu_val*180/math.pi)

elif mode == 2:
    nu_deg = float(input("nu (deg) = "))
    M = solve_M_from_nu(nu_deg*math.pi/180.0, gamma)
    print("Mach = ", M)

print("Done.")