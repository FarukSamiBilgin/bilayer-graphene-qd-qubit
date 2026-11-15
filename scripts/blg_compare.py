"""
Head-to-head dephasing comparison: single-dot operating point vs
double-well sweet spot. BOTH noise channels (detuning eps AND displacement
field D_ext) evaluated for the double well, so we don't hide an open channel.

Quasi-static inhomogeneous dephasing: sigma_Delta01 = sqrt( sum_k (sens_k * sigma_k)^2 )
T2* ratio = sigma_single / sigma_double.
"""
import numpy as np
from scipy.linalg import eigh
from blg_kp import Hkp, bulk_cond_min
from blg_dot import QX, QY, Npw
from blg_doubledot import confine_matrix, potential, dim
from blg_sweetspot import dot_split   # single dot (fast, verified)

# ---- single dot operating point ----
D_op, U0 = 0.167, 0.020
h = 0.004
s0_single = abs((dot_split(D_op+h,U0)-dot_split(D_op-h,U0))/(2*h))
d01_single = dot_split(D_op,U0)
print(f"SINGLE DOT  @ D_ext={D_op}, U0={1e3*U0}meV")
print(f"  Delta01={1e3*d01_single:.3f} meV | dDelta01/dD_ext = {s0_single:.4f} (linear)")

# ---- double well sweet spot ----
def Kin_of_D(D):
    K=np.zeros((dim,dim),dtype=complex)
    for a in range(Npw):
        K[4*a:4*a+4,4*a:4*a+4]=Hkp(QX[a],QY[a],D,+1)
    return K

_bulk={}
def gap(D):
    if D not in _bulk:_bulk[D]=bulk_cond_min(D)
    return _bulk[D]

def split_full(D, U0, Ldot, d, eps, mode, Cmat=None):
    if Cmat is None:
        Cmat = confine_matrix(potential(U0,Ldot,d,eps,mode))
    H = Kin_of_D(D) + Cmat
    E = eigh(H, eigvals_only=True)
    b=E[(E>1e-6)&(E<gap(D))]
    return (b[1]-b[0]) if len(b)>=2 else np.nan

Ldot, d, mode = 10.0, 30.0, 'depth'
heps = 0.02
# confinement at eps=0, +-heps (expensive: 3 builds)
C0  = confine_matrix(potential(U0,Ldot,d,0.0,mode))
Cp  = confine_matrix(potential(U0,Ldot,d,+heps,mode))
Cm  = confine_matrix(potential(U0,Ldot,d,-heps,mode))

d01_d = split_full(D_op,U0,Ldot,d,0.0,mode,C0)
# eps channel (at fixed D_op): 1st deriv ~0 by symmetry, 2nd deriv = curvature
sp = split_full(D_op,U0,Ldot,d,+heps,mode,Cp)
sm = split_full(D_op,U0,Ldot,d,-heps,mode,Cm)
d1_eps = (sp-sm)/(2*heps)
d2_eps = (sp-2*d01_d+sm)/heps**2
# D_ext channel (at eps=0, reuse C0): is it protected too?
sDp = split_full(D_op+h,U0,Ldot,d,0.0,mode,C0)
sDm = split_full(D_op-h,U0,Ldot,d,0.0,mode,C0)
s0_D_double = abs((sDp-sDm)/(2*h))

print(f"\nDOUBLE WELL @ d={d}nm, eps=0 (sweet spot), mode={mode}")
print(f"  Delta01={1e3*d01_d:.3f} meV")
print(f"  eps channel:   dDelta01/deps={d1_eps:.2e} (~0, symmetry) | "
      f"d2Delta01/deps2={d2_eps:.4f} eV (QUADRATIC)")
print(f"  D_ext channel: dDelta01/dD_ext={s0_D_double:.4f} "
      f"({'protected' if s0_D_double<0.3*s0_single else 'STILL OPEN'})")

# ---- dephasing comparison vs noise amplitude ----
# detuning eps is fractional depth asymmetry; physical detuning energy ~ U0*eps.
# Express eps-noise as fractional sigma_eps. D_ext-noise sigma_D in eV.
print("\n" + "="*64)
print("DEPHASING COMPARISON (T2* ratio = sigma_single / sigma_double)")
print("Assume charge noise drives BOTH eps (fractional) and D_ext (eV).")
print("="*64)
print(f"{'sigma_D(meV)':>12}{'sigma_eps':>10}{'sig_single':>11}{'sig_double':>11}{'T2*_gain':>10}")
for sigma_D_meV, sigma_eps in [(0.2,0.01),(0.2,0.05),(0.5,0.02),(0.5,0.05),(1.0,0.05),(1.0,0.10)]:
    sD = sigma_D_meV*1e-3
    sig_single = s0_single*sD
    # double: D_ext linear (open) + eps quadratic
    sig_double = np.sqrt((s0_D_double*sD)**2 + (0.5*abs(d2_eps)*sigma_eps**2)**2)
    print(f"{sigma_D_meV:>12.2f}{sigma_eps:>10.3f}"
          f"{1e3*sig_single:>11.4f}{1e3*sig_double:>11.4f}{sig_single/sig_double:>10.3f}")
