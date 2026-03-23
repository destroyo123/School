# Oblique Shock Calculator (Full Version)
# TI-84 CE Python compatible

import math

def deg2rad(d):
    return d*math.pi/180.0

def rad2deg(r):
    return r*180.0/math.pi

# Theta-Beta-Mach equation residual
def oblique_eq(beta, M1, delta, gamma):
    left = math.tan(delta)
    num = 2*(1/math.tan(beta))*(M1**2*math.sin(beta)**2 - 1)
    den = M1**2*(gamma + math.cos(2*beta)) + 2
    right = num/den
    return left - right

def solve_beta(M1, delta_deg, gamma, branch):

    delta = deg2rad(delta_deg)

    beta_min = math.asin(1.0/M1) + 1e-6
    beta_max = math.pi/2 - 1e-6

    roots = []
    N = 200
    step = (beta_max - beta_min)/N

    b1 = beta_min
    f1 = oblique_eq(b1, M1, delta, gamma)

    for i in range(N):
        b2 = b1 + step
        f2 = oblique_eq(b2, M1, delta, gamma)

        if f1*f2 < 0:
            low = b1
            high = b2
            for j in range(40):
                mid = 0.5*(low+high)
                fm = oblique_eq(mid, M1, delta, gamma)
                if oblique_eq(low, M1, delta, gamma)*fm < 0:
                    high = mid
                else:
                    low = mid
            roots.append(0.5*(low+high))

        b1 = b2
        f1 = f2

    if len(roots) == 0:
        return None

    if branch == 1:
        return min(roots)
    else:
        return max(roots)


# -------------------------
# USER INPUT
# -------------------------

gamma = float(input("Gamma (1.4 air) = "))
M1 = float(input("Upstream Mach M1 = "))
delta = float(input("Deflection angle delta (deg) = "))

print("1 = Weak shock")
print("2 = Strong shock")
branch = int(input("Branch #: "))

beta = solve_beta(M1, delta, gamma, branch)

if beta == None:
    print("No attached solution (shock detached).")
else:

    beta_deg = rad2deg(beta)

    # Normal Mach numbers
    M1n = M1*math.sin(beta)

    M2n = math.sqrt((1 + 0.5*(gamma-1)*M1n**2) /
                    (gamma*M1n**2 - 0.5*(gamma-1)))

    # Downstream Mach
    M2 = M2n / math.sin(beta - deg2rad(delta))

    # Static pressure ratio
    P2P1 = 1 + 2*gamma/(gamma+1)*(M1n**2 - 1)

    # Density ratio
    rho2rho1 = ((gamma+1)*M1n**2) / \
               ((gamma-1)*M1n**2 + 2)

    # Temperature ratio
    T2T1 = P2P1 / rho2rho1

    # Stagnation pressure ratio (normal shock relation)
    term1 = ((gamma+1)*M1n**2) / \
            ((gamma-1)*M1n**2 + 2)

    term2 = ((gamma+1)/(2*gamma*M1n**2 - (gamma-1)))

    P02P01 = term1**(gamma/(gamma-1)) * \
             term2**(1/(gamma-1))

    print("")
    print("===== Results =====")
    print("Beta (deg) = ", beta_deg)
    print("")
    print("M1n = ", M1n)
    print("M2n = ", M2n)
    print("M2 = ", M2)
    print("")
    print("rho2/rho1 = ", rho2rho1)
    print("P2/P1 = ", P2P1)
    print("T2/T1 = ", T2T1)
    print("P02/P01 = ", P02P01)

print("Done.")