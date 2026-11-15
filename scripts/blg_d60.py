"""
Fine scan around d=60nm: is the low D_ext-sensitivity a genuine protected
point or an orbital-crossing artifact? Track Delta01(d) AND lowest 3 bound
energies (to see if it's an anticrossing of the 1st excited state).
Also map D_ext-slope and eps-curvature across the window.
"""
import numpy as np
from scipy.linalg import eigh
from blg_kp import Hkp, bulk_cond_min
from blg_dot import QX, QY, Npw
from blg_doubledot import confine_matrix, potential, dim

U0, Ldot, mode = 0.020, 10.0, 'depth'
D_op, h = 0.167, 0.004
_bulk={}
def gap(D):
    if D not in _bulk: _bulk[D]=bulk_cond_min(D)
    return _bulk[D]
def Kin_of_D(D):
    K=np.zeros((dim,dim),dtype=complex)
    for a in range(Npw): K[4*a:4*a+4,4*a:4*a+4]=Hkp(QX[a],QY[a],D,+1)
    return K
def levels(D, Cmat, n=4):
    E=eigh(Kin_of_D(D)+Cmat, eigvals_only=True)
    b=E[(E>1e-6)&(E<gap(D))]
    return b[:n]
def split_full(D, Cmat):
    b=levels(D,Cmat,2)
    return (b[1]-b[0]) if len(b)>=2 else np.nan

print("Fine d-scan: bound-state structure + D_ext sensitivity")
print(f"{'d(nm)':>6}{'E0':>8}{'E1':>8}{'E2':>8}{'D01(meV)':>9}{'D12(meV)':>9}{'|dD01/dDext|':>13}")
for d in np.arange(52.,69.,2.0):
    C0=confine_matrix(potential(U0,Ldot,d,0.0,mode))
    b=levels(D_op,C0,3)
    s0=b[1]-b[0]
    s12=(b[2]-b[1]) if len(b)>=3 else np.nan
    sl=abs((split_full(D_op+h,C0)-split_full(D_op-h,C0))/(2*h))
    print(f"{d:>6.0f}{1e3*b[0]:>8.2f}{1e3*b[1]:>8.2f}"
          f"{1e3*b[2] if len(b)>=3 else float('nan'):>8.2f}"
          f"{1e3*s0:>9.3f}{1e3*s12:>9.3f}{sl:>13.4f}")
