"""Fast AC test on the real BLG dot (optimized)."""
import numpy as np
from scipy.linalg import eigh
from blg_kp import Hkp
from blg_dot import QX, QY, Npw, Vft, L_box

U0, Ldot = 0.020, 15.0
D_op = 0.167

# --- Precompute D-independent pieces once ---
dim = 4*Npw
# kinetic at D=0
H0 = np.zeros((dim, dim), dtype=complex)
for a in range(Npw):
    H0[4*a:4*a+4, 4*a:4*a+4] = Hkp(QX[a], QY[a], 0.0, +1)
# confinement (D-independent), scalar x I4
for a in range(Npw):
    for b in range(Npw):
        val = Vft(QX[a]-QX[b], QY[a]-QY[b], U0, Ldot)
        if val != 0.0:
            for o in range(4):
                H0[4*a+o, 4*b+o] += val
# D-operator: +1/2 on A1,B1 (o=0,1); -1/2 on A2,B2 (o=2,3)
Ddiag = np.zeros(dim)
for a in range(Npw):
    Ddiag[4*a+0] = +0.5; Ddiag[4*a+1] = +0.5
    Ddiag[4*a+2] = -0.5; Ddiag[4*a+3] = -0.5

_cache = {}
def d01(D):
    key = round(D, 7)
    if key in _cache: return _cache[key]
    H = H0 + np.diag(D*Ddiag).astype(complex)
    # in-gap positive states only
    ev = eigh(H, eigvals_only=True, subset_by_value=[1e-5, 0.12])
    val = ev[1]-ev[0]
    _cache[key] = val
    return val

NPH = 12
ph = np.linspace(0, 2*np.pi, NPH, endpoint=False)
cosph = np.cos(ph)

def periodavg(D, A):
    if A == 0.0: return d01(D)
    return np.mean([d01(D + A*c) for c in cosph])

def slope(A, h=0.004):
    return (periodavg(D_op+h, A) - periodavg(D_op-h, A)) / (2*h)

if __name__ == "__main__":
    import sys
    A_list = [float(x) for x in sys.argv[1:]] or [0.0, 0.030]
    print(f"dim={dim}, NPH={NPH}, D_op={D_op}")
    print(f"static Delta01 = {1e3*d01(D_op):.4f} meV")
    s0 = abs(slope(0.0))
    print(f"{'A(meV)':>7} {'<D01>_T':>9} {'|slope|':>9} {'sig_ratio':>10} {'T2*_ratio':>10}")
    for A in A_list:
        avg = periodavg(D_op, A); s = abs(slope(A))
        print(f"{1e3*A:>7.1f} {1e3*avg:>9.4f} {s:>9.4f} {s/s0:>10.4f} {s0/s:>10.4f}")
    print(f"anchor |avg(0)-static| = {abs(periodavg(D_op,0.0)-d01(D_op)):.2e}")
