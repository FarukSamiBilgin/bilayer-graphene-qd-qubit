"""Test: is dDelta01/dD_ext irreducible in BLG dots? Scan double-well d and
single-dot, look for ANY config with small D_ext-slope at useful splitting."""
import numpy as np
from scipy.linalg import eigh
from blg_kp import Hkp, bulk_cond_min
from blg_dot import QX, QY, Npw
from blg_doubledot import confine_matrix, potential, dim

U0, Ldot, mode = 0.020, 10.0, 'depth'
D_op = 0.167; h=0.004
_bulk={}
def gap(D):
    if D not in _bulk:_bulk[D]=bulk_cond_min(D)
    return _bulk[D]
def Kin_of_D(D):
    K=np.zeros((dim,dim),dtype=complex)
    for a in range(Npw): K[4*a:4*a+4,4*a:4*a+4]=Hkp(QX[a],QY[a],D,+1)
    return K
def split_full(D, Cmat):
    E=eigh(Kin_of_D(D)+Cmat, eigvals_only=True)
    b=E[(E>1e-6)&(E<gap(D))]
    return (b[1]-b[0]) if len(b)>=2 else np.nan

print("Double-well: D_ext sensitivity at eps=0 sweet spot, vs separation d")
print(f"{'d(nm)':>6}{'Delta01(meV)':>13}{'|dD01/dD_ext|':>15}{'ratio_to_gap':>13}")
for d in (20.,25.,30.,40.,50.,60.,80.):
    C0=confine_matrix(potential(U0,Ldot,d,0.0,mode))
    s0=split_full(D_op,C0)
    sl=abs((split_full(D_op+h,C0)-split_full(D_op-h,C0))/(2*h))
    # dimensionless sensitivity relative to splitting: (dD01/dD_ext)/(D01) -- fractional
    frac = sl/s0 if s0>0 else np.nan
    print(f"{d:>6.0f}{1e3*s0:>13.3f}{sl:>15.4f}{frac:>13.2f}")
