"""
Double-well BLG dot in the plane-wave basis. Two Gaussians at x=+-d/2.
Two detuning definitions:
  mode='depth' : eps tilts well DEPTH  (left -U0(1+e/2), right -U0(1-e/2))
  mode='tilt'  : add linear potential eps * x / x_scale  (lateral field)
Shift theorem: FT of Gaussian at x0 gains phase exp(-i qx x0).

Sweet-spot anchor: at eps=0 the two wells are identical -> by symmetry
dDelta01/deps MUST vanish (this is the candidate sweet spot). If it does,
the second derivative sets the residual sensitivity & whether it's useful.
"""
import numpy as np
from scipy.linalg import eigh
from blg_kp import Hkp, bulk_cond_min
from blg_dot import QX, QY, Npw, L_box

dim = 4*Npw
D_FIX = 0.167          # fixed displacement-field gap (operating point)

# kinetic at D=D_FIX (fixed throughout the eps scan)
Kin = np.zeros((dim, dim), dtype=complex)
for a in range(Npw):
    Kin[4*a:4*a+4, 4*a:4*a+4] = Hkp(QX[a], QY[a], D_FIX, +1)

# linear-x operator in plane-wave basis: <q|x|q'> -> off-diagonal; easier to
# build x as a dense matrix via FFT-free approach is messy. Instead implement
# 'tilt' as a smooth odd potential V_tilt(r)=tanh(x/w): its FT we sample numerically.
# Simpler & robust: build confinement DIRECTLY in real space on a grid, then
# project onto plane waves. Use a real-space grid and FFT.

# real-space grid
Ng = 64
xs = (np.arange(Ng)-Ng/2)*(L_box/Ng)
X, Y = np.meshgrid(xs, xs, indexing='ij')
# plane-wave indices -> map each (qx,qy) to integer (mx,my)
dq = 2*np.pi/L_box
MX = np.round(QX/dq).astype(int) % Ng
MY = np.round(QY/dq).astype(int) % Ng

def confine_matrix(Vr):
    """Project a real-space scalar potential Vr(X,Y) onto the plane-wave basis."""
    Vk = np.fft.fft2(Vr)/(Ng*Ng)          # Vk[m] = (1/N) sum V e^{-i 2pi m n/N}
    C = np.zeros((dim, dim), dtype=complex)
    for a in range(Npw):
        for b in range(Npw):
            mx = (MX[a]-MX[b]) % Ng
            my = (MY[a]-MY[b]) % Ng
            val = Vk[mx, my]
            if val != 0:
                for o in range(4):
                    C[4*a+o, 4*b+o] += val
    return C

def potential(U0, Ldot, d, eps, mode):
    xl, xr = -d/2, +d/2
    gl = np.exp(-((X-xl)**2+Y**2)/(2*Ldot**2))
    gr = np.exp(-((X-xr)**2+Y**2)/(2*Ldot**2))
    if mode=='depth':
        return -U0*((1+eps/2)*gl + (1-eps/2)*gr)
    elif mode=='tilt':
        return -U0*(gl+gr) + eps*np.tanh(X/(d/2+1e-9))*U0  # odd tilt ~ +-U0*eps at wells
    else:
        raise ValueError

_bulk = {}
def gap():
    if D_FIX not in _bulk: _bulk[D_FIX]=bulk_cond_min(D_FIX)
    return _bulk[D_FIX]

def split(U0, Ldot, d, eps, mode):
    H = Kin + confine_matrix(potential(U0,Ldot,d,eps,mode))
    E = eigh(H, eigvals_only=True)
    cmin = gap()
    bound = E[(E>1e-6)&(E<cmin)]
    return (bound[1]-bound[0]) if len(bound)>=2 else np.nan

if __name__ == "__main__":
    U0, Ldot = 0.020, 10.0
    print(f"Double-well BLG dot, D_ext={D_FIX}, U0={1e3*U0}meV, Ldot={Ldot}nm, Ng={Ng}")
    for mode in ('depth','tilt'):
        print(f"\n=== mode={mode} ===")
        for d in (30.0, 50.0, 80.0):
            s0 = split(U0,Ldot,d,0.0,mode)
            he = 0.02
            sp = split(U0,Ldot,d,+he,mode); sm = split(U0,Ldot,d,-he,mode)
            d1 = (sp-sm)/(2*he)                 # dDelta01/deps at eps=0 (ANCHOR: ~0)
            d2 = (sp-2*s0+sm)/he**2
            print(f" d={d:4.0f}nm: Delta01(eps=0)={1e3*s0:7.3f} meV | "
                  f"dD/deps={d1: .2e} (anchor~0) | d2D/deps2={d2: .3f}")
