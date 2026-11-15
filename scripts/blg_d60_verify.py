"""At d=60nm: verify eps sweet spot AND fair head-to-head vs single dot,
using quality factor Q ~ Delta01/sigma_Delta01 (penalizes small splitting)."""
import numpy as np
from scipy.linalg import eigh
from blg_kp import Hkp, bulk_cond_min
from blg_dot import QX, QY, Npw
from blg_doubledot import confine_matrix, potential, dim
from blg_sweetspot import dot_split

U0, Ldot, mode, d = 0.020, 10.0, 'depth', 60.0
D_op, h, heps = 0.167, 0.004, 0.02
_bulk={}
def gap(D):
    if D not in _bulk:_bulk[D]=bulk_cond_min(D)
    return _bulk[D]
def Kin(D):
    K=np.zeros((dim,dim),dtype=complex)
    for a in range(Npw): K[4*a:4*a+4,4*a:4*a+4]=Hkp(QX[a],QY[a],D,+1)
    return K
def spl(D,C):
    b=eigh(Kin(D)+C,eigvals_only=True); b=b[(b>1e-6)&(b<gap(D))]
    return (b[1]-b[0]) if len(b)>=2 else np.nan

C0=confine_matrix(potential(U0,Ldot,d,0.0,mode))
Cp=confine_matrix(potential(U0,Ldot,d,+heps,mode))
Cm=confine_matrix(potential(U0,Ldot,d,-heps,mode))
d01=spl(D_op,C0)
d1e=(spl(D_op,Cp)-spl(D_op,Cm))/(2*heps)
d2e=(spl(D_op,Cp)-2*d01+spl(D_op,Cm))/heps**2
sD =abs((spl(D_op+h,C0)-spl(D_op-h,C0))/(2*h))

# single dot reference
d01_s=dot_split(D_op,U0); sD_s=abs((dot_split(D_op+h,U0)-dot_split(D_op-h,U0))/(2*h))

print(f"DOUBLE WELL d={d}nm, eps=0:")
print(f"  Delta01={1e3*d01:.3f} meV")
print(f"  eps: dD/deps={d1e:.2e} (~0 sweet spot), d2D/deps2={d2e:.4f} eV")
print(f"  D_ext: |dD/dDext|={sD:.4f}")
print(f"\nSINGLE DOT D_ext={D_op}: Delta01={1e3*d01_s:.3f} meV, |dD/dDext|={sD_s:.4f}")
print("\n"+"="*58)
print("Realistic charge noise: sigma_eps small (mu-eV detuning), D_ext dominant")
print("Quality factor Q = Delta01 / sigma_Delta01  (ops before dephasing)")
print("="*58)
print(f"{'sigma_D(meV)':>12}{'Q_single':>10}{'Q_double':>10}{'Q_gain':>8}")
for sD_meV in (0.2,0.5,1.0):
    s=sD_meV*1e-3
    sig_s=sD_s*s
    sig_d=np.sqrt((sD*s)**2+(0.5*abs(d2e)*0.002**2)**2)  # eps quadratic, sigma_eps=0.002
    Qs=d01_s/sig_s; Qd=d01/sig_d
    print(f"{sD_meV:>12.2f}{Qs:>10.1f}{Qd:>10.1f}{Qd/Qs:>8.2f}")
