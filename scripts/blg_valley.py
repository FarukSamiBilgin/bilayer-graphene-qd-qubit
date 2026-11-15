"""
NEW DIRECTION: valley qubit charge-noise sensitivity.
Orbital magnetic moment mu_orb of the dot ground state, per valley.
Valley splitting in small B: Delta_v = 2 mu_orb B (valleys have opposite mu).
Hypothesis: d(mu_orb)/dD_ext is small -> valley qubit is charge-noise protected,
unlike the orbital qubit (fractional slope ~28/eV).

mu_orb ~ <psi| (x v_y - y v_x) |psi>, v_i = dH/dq_i (band-space matrices).
Computed in real space via FFT of plane-wave coefficients (no q-grid artifacts).
ANCHOR: mu_orb(K) = -mu_orb(K')  (time-reversal symmetry).
"""
import numpy as np
from scipy.linalg import eigh
from blg_kp import Hkp, bulk_cond_min, v, v3
from blg_dot import QX, QY, Npw, Vft, L_box

Ldot = 15.0
dim = 4*Npw

# velocity matrices Vx=dH/dqx, Vy=dH/dqy (band space, depend on valley xi)
def velocity(xi):
    Vx = np.array([[0,    v*xi, 0,    v3*xi],
                   [v*xi, 0,    0,    0   ],
                   [0,    0,    0,    v*xi],
                   [v3*xi,0,    v*xi, 0   ]], dtype=complex)
    Vy = np.array([[0,    -1j*v, 0,    1j*v3],
                   [1j*v, 0,     0,    0    ],
                   [0,    0,     0,    -1j*v],
                   [-1j*v3,0,    1j*v, 0    ]], dtype=complex)
    return Vx, Vy

# real-space grid for wavefunction reconstruction
Ng = 96
xs = (np.arange(Ng)-Ng/2)*(L_box/Ng)
Xg, Yg = np.meshgrid(xs, xs, indexing='ij')
dq = 2*np.pi/L_box
MX = np.round(QX/dq).astype(int) % Ng
MY = np.round(QY/dq).astype(int) % Ng

def build_H(D, U0, xi):
    H = np.zeros((dim,dim),dtype=complex)
    for a in range(Npw):
        H[4*a:4*a+4,4*a:4*a+4]=Hkp(QX[a],QY[a],D,xi)
    for a in range(Npw):
        for b in range(Npw):
            val=Vft(QX[a]-QX[b],QY[a]-QY[b],U0,Ldot)
            if val!=0:
                for o in range(4): H[4*a+o,4*b+o]+=val
    return H

def ground_state(D,U0,xi):
    H=build_H(D,U0,xi)
    E,W=eigh(H)
    cmin=bulk_cond_min(D)
    idx=np.where((E>1e-6)&(E<cmin))[0]
    return W[:,idx[0]]   # ground bound state, 4*Npw coeffs

def psi_realspace(c):
    """Reconstruct 4 band components on the Ng grid from plane-wave coeffs."""
    psi=np.zeros((4,Ng,Ng),dtype=complex)
    for o in range(4):
        F=np.zeros((Ng,Ng),dtype=complex)
        for a in range(Npw):
            F[MX[a],MY[a]] += c[4*a+o]
        psi[o]=np.fft.ifft2(F)*(Ng*Ng)   # to real space
    # normalize on grid
    nrm=np.sqrt(np.sum(np.abs(psi)**2))
    return psi/nrm

def mu_orb(D,U0,xi):
    c=ground_state(D,U0,xi)
    psi=psi_realspace(c)               # (4,Ng,Ng)
    Vx,Vy=velocity(xi)
    # local currents J_i(r)=Re[ psi^dag V_i psi ] at each r
    P=psi.reshape(4,-1)                # 4 x Npts
    Jx=np.real(np.einsum('ir,ij,jr->r',P.conj(),Vx,P))
    Jy=np.real(np.einsum('ir,ij,jr->r',P.conj(),Vy,P))
    x=Xg.ravel(); y=Yg.ravel()
    m = np.sum(x*Jy - y*Jx)            # reduced orbital moment (drop e/2hbar)
    return m

if __name__ == "__main__":
    U0=0.020
    print("Orbital magnetic moment mu_orb of dot ground state (reduced units)")
    print(f"{'D_ext':>7}{'mu(K)':>12}{'mu(Kp)':>12}{'sum(anchor~0)':>15}")
    for D in (0.10,0.15,0.167,0.20):
        mK=mu_orb(D,U0,+1); mKp=mu_orb(D,U0,-1)
        print(f"{D:>7.3f}{mK:>12.4f}{mKp:>12.4f}{mK+mKp:>15.2e}")

if __name__ == "__main__":
    print("\n" + "="*60)
    print("VALLEY vs ORBITAL qubit: fractional D_ext sensitivity")
    print("="*60)
    U0=0.020; D=0.167; h=0.004
    # valley: fractional d(mu)/dD / mu  (since Delta_v = 2 mu B)
    mu0=mu_orb(D,U0,+1)
    dmu=(mu_orb(D+h,U0,+1)-mu_orb(D-h,U0,+1))/(2*h)
    frac_valley=abs(dmu/mu0)
    # orbital: from earlier, Delta01 and slope
    from blg_sweetspot import dot_split
    d01=dot_split(D,U0); s01=abs((dot_split(D+h,U0)-dot_split(D-h,U0))/(2*h))
    frac_orbital=s01/d01
    print(f"  ORBITAL qubit: Delta01={1e3*d01:.3f} meV, "
          f"frac sens = {frac_orbital:.2f} /eV")
    print(f"  VALLEY qubit:  mu_orb={mu0:.3f}, "
          f"frac sens = {frac_valley:.2f} /eV")
    print(f"  => VALLEY is {frac_orbital/frac_valley:.2f}x LESS charge-sensitive")

    # convergence of mu_orb fractional sensitivity vs grid Ng
    print("\nConvergence of valley fractional sensitivity vs Ng grid:")
    import blg_valley as BV
    for ng in (64,96,128):
        BV.Ng=ng
        BV.xs=(np.arange(ng)-ng/2)*(L_box/ng)
        BV.Xg,BV.Yg=np.meshgrid(BV.xs,BV.xs,indexing='ij')
        BV.MX=np.round(QX/dq).astype(int)%ng
        BV.MY=np.round(QY/dq).astype(int)%ng
        m0=BV.mu_orb(D,U0,+1); dm=(BV.mu_orb(D+h,U0,+1)-BV.mu_orb(D-h,U0,+1))/(2*h)
        print(f"   Ng={ng}: mu={m0:.4f}, frac_sens={abs(dm/m0):.3f} /eV")
