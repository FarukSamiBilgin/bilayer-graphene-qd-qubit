"""
Confined BLG quantum dot: 4-band k.p in a plane-wave basis with a
Gaussian scalar well V(r) = -U0 exp(-r^2 / (2 L_dot^2)).
Reproduce Table 1 (Delta_ext sweep -> Delta01) and test curvature.
"""
import numpy as np
from blg_kp import Hkp, g1, bulk_cond_min

# ---- plane-wave box ----
L_box = 120.0                  # nm
dq    = 2*np.pi/L_box          # nm^-1
qmax  = 0.5                    # nm^-1 cutoff

# build plane-wave momentum list |q| < qmax
ns = np.arange(-int(qmax/dq)-1, int(qmax/dq)+2)
QX, QY = [], []
for i in ns:
    for j in ns:
        qx, qy = i*dq, j*dq
        if qx*qx + qy*qy < qmax*qmax:
            QX.append(qx); QY.append(qy)
QX = np.array(QX); QY = np.array(QY)
Npw = len(QX)
print(f"L_box={L_box} nm, dq={dq:.4f}, qmax={qmax} -> N_pw={Npw}, dim={4*Npw}")

def Vft(qx, qy, U0, Ldot):
    # 2D FT of -U0 exp(-r^2/(2 Ldot^2)), normalized by box area
    return -U0 * 2*np.pi*Ldot**2 * np.exp(-(qx**2+qy**2)*Ldot**2/2) / L_box**2

def build_H(D, U0, Ldot, xi=+1):
    dim = 4*Npw
    H = np.zeros((dim, dim), dtype=complex)
    # kinetic: block-diagonal k.p at each plane wave
    for a in range(Npw):
        H[4*a:4*a+4, 4*a:4*a+4] = Hkp(QX[a], QY[a], D, xi)
    # confinement: scalar V_FT(q-q') x Identity_4
    for a in range(Npw):
        for b in range(Npw):
            val = Vft(QX[a]-QX[b], QY[a]-QY[b], U0, Ldot)
            if val != 0.0:
                for o in range(4):
                    H[4*a+o, 4*b+o] += val
    return H

def dot_states(D, U0, Ldot, xi=+1):
    H = build_H(D, U0, Ldot, xi)
    E = np.linalg.eigvalsh(H)
    cmin = bulk_cond_min(D)            # bulk conduction edge
    # bound states: in-gap, below bulk conduction min, above 0
    bound = E[(E > 0) & (E < cmin)]
    return bound, cmin

if __name__ == "__main__":
    U0, Ldot = 0.020, 15.0     # 20 meV, 15 nm
    print(f"\nTABLE 1 reproduction: U0={1e3*U0} meV, L_dot={Ldot} nm")
    print(f"{'D_ext':>7} {'bulk_gap':>9} {'E0(meV)':>9} {'D01(meV)':>9} {'#bound':>7}")
    print("-"*48)
    rows = []
    for D in (0.05, 0.10, 0.15, 0.167, 0.20, 0.25, 0.30):
        bound, cmin = dot_states(D, U0, Ldot)
        if len(bound) >= 2:
            E0, E1 = bound[0], bound[1]
            d01 = E1 - E0
            rows.append((D, d01))
            print(f"{D:>7.3f} {1e3*cmin:>9.2f} {1e3*E0:>9.2f} "
                  f"{1e3*d01:>9.3f} {len(bound):>7d}")
        else:
            print(f"{D:>7.3f} {1e3*cmin:>9.2f} {'--':>9} {'--':>9} {len(bound):>7d}")

    # curvature test at operating point
    print("\nCURVATURE of Delta01(D_ext) around D_ext=0.167:")
    Ds = np.array([0.150, 0.167, 0.184])
    d01s = []
    for D in Ds:
        b,_ = dot_states(D, U0, Ldot)
        d01s.append(b[1]-b[0])
    d01s = np.array(d01s)
    first  = (d01s[2]-d01s[0])/(Ds[2]-Ds[0])
    second = (d01s[2]-2*d01s[1]+d01s[0])/((Ds[1]-Ds[0])**2)
    print(f"  Delta01 at {Ds} = {[f'{1e3*x:.3f}' for x in d01s]} meV")
    print(f"  d(Delta01)/dD_ext   ~ {first:.4f}  (sign of slope)")
    print(f"  d2(Delta01)/dD_ext2 ~ {second:.4f}  "
          f"({'CONCAVE' if second<0 else 'CONVEX'})")
