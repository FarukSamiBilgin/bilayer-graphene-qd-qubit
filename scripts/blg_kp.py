"""
4-band k.p (McCann-Koshino) bilayer graphene solver.
Basis: (A1, B1, A2, B2).  A1,B1 = bottom layer (+D/2); A2,B2 = top (-D/2).
gamma1 couples dimer sites B1-A2; gamma3 (trigonal) couples A1-B2.
Units: energy eV, length nm, momentum nm^-1.
"""
import numpy as np

# ---- parameters (paper values) ----
a   = 0.246      # nm
g0  = 3.16       # eV  (intralayer)
g1  = 0.38       # eV  (interlayer dimer)
g3  = 0.38       # eV  (trigonal)
v   = (np.sqrt(3)/2) * a * g0     # eV*nm  ~0.673
v3  = (np.sqrt(3)/2) * a * g3     # eV*nm  ~0.081

def Hkp(qx, qy, D, xi=+1):
    """4x4 k.p Hamiltonian at momentum (qx,qy), gap D, valley xi=+-1."""
    pi  = xi*qx + 1j*qy
    pid = xi*qx - 1j*qy
    H = np.array([
        [ +D/2,      v*pid,     0.0,       v3*pi ],
        [ v*pi,      +D/2,      g1,        0.0   ],
        [ 0.0,       g1,        -D/2,      v*pid ],
        [ v3*pid,    0.0,       v*pi,      -D/2  ],
    ], dtype=complex)
    return H

# ============ ANCHOR 1: Hermiticity ============
def anchor1():
    worst = 0.0
    rng = np.random.default_rng(0)
    for _ in range(200):
        qx,qy = rng.uniform(-0.5,0.5,2)
        D = rng.uniform(0,0.3)
        for xi in (+1,-1):
            H = Hkp(qx,qy,D,xi)
            worst = max(worst, np.abs(H - H.conj().T).max())
    return worst

# ============ ANCHOR 2: bulk conduction-band minimum ============
def bulk_cond_min(D, xi=+1, qmax=0.15, n=401):
    qs = np.linspace(-qmax, qmax, n)
    cmin = np.inf
    for qx in qs:
        for qy in qs:
            E = np.linalg.eigvalsh(Hkp(qx,qy,D,xi))
            # lowest positive (conduction) eigenvalue
            pos = E[E > 0]
            if len(pos): cmin = min(cmin, pos.min())
    return cmin

def anchor2():
    out = []
    for D in (0.05, 0.10, 0.167, 0.30):
        num = bulk_cond_min(D)
        ana = g1*D/(2*np.sqrt(g1**2 + D**2))   # gamma3=0 analytic
        out.append((D, num, ana, 100*(num-ana)/ana))
    return out

# ============ ANCHOR 3: C3 symmetry of trigonal warping ============
def anchor3(D=0.10):
    # sample conduction energy on a ring at fixed |q|, check 3-fold periodicity
    qr = 0.08
    th = np.linspace(0, 2*np.pi, 360, endpoint=False)
    E = np.array([np.linalg.eigvalsh(Hkp(qr*np.cos(t), qr*np.sin(t), D))[2]
                  for t in th])
    # FFT: with C3 warping, dominant angular harmonic is m=3
    amp = np.abs(np.fft.rfft(E - E.mean()))
    m_dom = int(np.argmax(amp))
    # also compare g3=0 (should be ~isotropic) by temporarily zeroing v3
    return m_dom, amp[m_dom], E.std()

def anchor3_iso():
    global v3
    save = v3; v3 = 0.0
    qr=0.08; th=np.linspace(0,2*np.pi,360,endpoint=False)
    E = np.array([np.linalg.eigvalsh(Hkp(qr*np.cos(t),qr*np.sin(t),0.10))[2] for t in th])
    v3 = save
    return E.std()

# ============ ANCHOR 4: D=0 gapless ============
def anchor4():
    # gap at/near q=0 with D=0 should vanish
    gaps = []
    for q in (0.0, 0.001, 0.01):
        E = np.linalg.eigvalsh(Hkp(q,0.0,0.0))
        gaps.append(E[2]-E[1])
    return gaps

# ============ ANCHOR 5: K vs K' identical spectra ============
def anchor5(D=0.167):
    rng = np.random.default_rng(1)
    worst = 0.0
    for _ in range(300):
        qx,qy = rng.uniform(-0.4,0.4,2)
        Ep = np.linalg.eigvalsh(Hkp(qx,qy,D,+1))
        # TRS maps K(q) -> K'(-q)
        Em = np.linalg.eigvalsh(Hkp(-qx,-qy,D,-1))
        worst = max(worst, np.abs(np.sort(Ep)-np.sort(Em)).max())
    return worst

if __name__ == "__main__":
    print("v  =", round(v,4), "eV*nm   v3 =", round(v3,4), "eV*nm")
    print("="*60)
    print("ANCHOR 1 (Hermiticity): max|H-H^dag| =", f"{anchor1():.2e}")
    print("="*60)
    print("ANCHOR 2 (bulk conduction min vs analytic gamma3=0):")
    print(f"{'D (eV)':>8} {'numeric(meV)':>13} {'analytic(meV)':>14} {'diff %':>8}")
    for D,num,ana,pct in anchor2():
        print(f"{D:>8.3f} {1e3*num:>13.2f} {1e3*ana:>14.2f} {pct:>8.1f}")
    print("="*60)
    m_dom, amp, std = anchor3()
    print(f"ANCHOR 3 (C3 trigonal): dominant angular harmonic m={m_dom} "
          f"(expect 3); ring energy std={1e3*std:.3f} meV")
    print(f"          gamma3=0 control: ring std={1e3*anchor3_iso():.4f} meV (expect ~0)")
    print("="*60)
    print("ANCHOR 4 (D=0 gapless): gaps near q=0 =",
          [f"{1e3*g:.3f}" for g in anchor4()], "meV")
    print("="*60)
    print(f"ANCHOR 5 (K vs K' Kramers): max spectral diff = {anchor5():.2e}")
