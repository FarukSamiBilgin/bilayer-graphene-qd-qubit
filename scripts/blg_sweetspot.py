"""
Static 2D sweet-spot search for the BLG valley qubit.
Charge noise couples to D_ext. A sweet spot = point where
d(Delta01)/d(D_ext) = 0 at a USEFUL (non-vanishing) splitting.
Second control axis: dot depth U0 (separate gate). Scan (D_ext,U0).

Speed: confinement ~ U0 (linear), kinetic split into D=0 part + D*Ddiag.
So H(D,U0) = Kin0 + D*Dmat + U0*Cmat, precomputed once; each grid pt = 1 eigh.
"""
import numpy as np
from scipy.linalg import eigh
from blg_kp import Hkp, bulk_cond_min
from blg_dot import QX, QY, Npw, Vft, L_box

Ldot = 15.0
dim = 4*Npw

# --- precompute the three building blocks ONCE ---
Kin0 = np.zeros((dim, dim), dtype=complex)
for a in range(Npw):
    Kin0[4*a:4*a+4, 4*a:4*a+4] = Hkp(QX[a], QY[a], 0.0, +1)

Dmat = np.zeros(dim)                      # layer-asymmetry diagonal (the D_ext operator)
for a in range(Npw):
    Dmat[4*a+0]=+0.5; Dmat[4*a+1]=+0.5; Dmat[4*a+2]=-0.5; Dmat[4*a+3]=-0.5
Dmat = np.diag(Dmat).astype(complex)

Cmat = np.zeros((dim, dim), dtype=complex)  # confinement at U0=1 (scales linearly)
for a in range(Npw):
    for b in range(Npw):
        val = Vft(QX[a]-QX[b], QY[a]-QY[b], 1.0, Ldot)   # U0=1
        if val != 0.0:
            for o in range(4):
                Cmat[4*a+o, 4*b+o] += val

_bulkcache = {}
def gap(D):
    k=round(D,6)
    if k not in _bulkcache: _bulkcache[k]=bulk_cond_min(D)
    return _bulkcache[k]

def dot_split(D, U0):
    H = Kin0 + D*Dmat + U0*Cmat
    E = eigh(H, eigvals_only=True)
    cmin = gap(D)
    bound = E[(E>1e-6) & (E<cmin)]
    return (bound[1]-bound[0]) if len(bound)>=2 else np.nan

if __name__ == "__main__":
    # VERIFY: U0=0.020 slice must reproduce Table 1
    print("VERIFY vs Table 1 (U0=20 meV):")
    for D in (0.05,0.10,0.15,0.167,0.20):
        print(f"  D={D:.3f}  Delta01={1e3*dot_split(D,0.020):.3f} meV")
    print("(expect 5.130, 3.471, 1.524, 0.989, 0.347)\n")

    # 2D scan: find d(Delta01)/dD_ext = 0 at useful splitting
    Ds = np.linspace(0.08, 0.22, 15)
    U0s = np.array([0.012,0.016,0.020,0.024,0.028,0.032,0.040])
    h = 0.004
    print("2D sweet-spot scan: |d Delta01/d D_ext| (meV/eV), * = local min in D at useful gap")
    print("U0\\D   " + "  ".join(f"{1e3*D:5.0f}" for D in Ds))
    best = []
    for U0 in U0s:
        row=[]
        for D in Ds:
            sp = dot_split(D,U0)
            dp = (dot_split(D+h,U0)-dot_split(D-h,U0))/(2*h)
            row.append((D,sp,abs(dp)))
        # find min |slope| among points with useful splitting (>0.3 meV)
        useful=[(D,sp,ds) for (D,sp,ds) in row if sp>3e-4]
        line=f"{1e3*U0:4.0f}  "
        for (D,sp,ds) in row:
            mark="*" if (sp>3e-4 and ds==min(x[2] for x in useful)) else " "
            line+=f"{ds:6.3f}{mark}"
        print(line)
        if useful:
            bm=min(useful,key=lambda x:x[2])
            best.append((U0,)+bm)
    print("\nBest (lowest |slope|) per U0, at useful splitting:")
    print(f"{'U0(meV)':>8}{'D_ext':>8}{'D01(meV)':>10}{'|slope|':>9}")
    for U0,D,sp,ds in best:
        print(f"{1e3*U0:>8.0f}{D:>8.3f}{1e3*sp:>10.3f}{ds:>9.3f}")
# (Extracted from notebooks/birlesimblg.ipynb, cell 4.)
