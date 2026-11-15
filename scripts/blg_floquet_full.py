"""
Full-Hilbert Floquet test of the period-averaged AC-protection claim.

Evolves the two lowest bound (qubit) states of the BLG dot through one drive
period of D_ext(t) = D_op + A cos(omega t) in the complete 1172-dim plane-wave
Hilbert space (Strang/Trotter splitting in the static eigenbasis, no
projection) and reports (i) leakage out of the qubit subspace and (ii) the
Floquet quasi-energy splitting.  omega = 6.6 * Delta01 (fast-drive regime
assumed by the period-averaging method).

Anchors: A = 0 -> leakage ~ 0 and Floquet splitting == static Delta01;
time-step convergence (nstep = 200/300/500).
(Extracted from the notebook notebooks/birlesimblg.ipynb, cell 3.)
"""
import numpy as np
from scipy.linalg import eigh
from blg_ac_fast import H0, Ddiag, d01

D_op = 0.167

def setup(D_c):
    H = H0 + np.diag(D_c*Ddiag).astype(complex)
    E, W = eigh(H)
    Vt = (W.conj().T * Ddiag) @ W
    mu, Q = eigh(Vt)
    Qd = Q.conj().T
    pos = np.where(E > 1e-6)[0]
    return E, mu, Q, Qd, pos[0], pos[1]

def evolve_period(E, mu, Q, Qd, psi0, A, omega, nstep):
    T = 2*np.pi/omega; dt = T/nstep
    free = np.exp(-1j*E*dt/2)
    psi = psi0.copy()
    for k in range(nstep):
        psi = free*psi
        c = A*np.cos(omega*(k+0.5)*dt)
        psi = Q @ (np.exp(-1j*dt*c*mu) * (Qd @ psi))
        psi = free*psi
    return psi, T

def floquet_full(D_c, A, omega, nstep=300):
    E, mu, Q, Qd, q0, q1 = setup(D_c)
    dim = len(E)
    e0 = np.zeros(dim, dtype=complex); e0[q0]=1
    e1 = np.zeros(dim, dtype=complex); e1[q1]=1
    p0, T = evolve_period(E, mu, Q, Qd, e0, A, omega, nstep)
    p1, _ = evolve_period(E, mu, Q, Qd, e1, A, omega, nstep)
    M = np.array([[np.vdot(e0,p0), np.vdot(e0,p1)],
                  [np.vdot(e1,p0), np.vdot(e1,p1)]])
    leak0 = 1-(abs(np.vdot(e0,p0))**2+abs(np.vdot(e1,p0))**2)
    leak1 = 1-(abs(np.vdot(e0,p1))**2+abs(np.vdot(e1,p1))**2)
    leak = max(leak0, leak1)
    Uq,s,Vh = np.linalg.svd(M); Mu = Uq@Vh
    ph = np.angle(np.linalg.eigvals(Mu))
    dphi = abs(((ph[0]-ph[1]+np.pi)%(2*np.pi))-np.pi)
    return leak, dphi/T

if __name__ == "__main__":
    d0 = d01(D_op); omega = 6.6*d0
    print(f"static Delta01={1e3*d0:.4f} meV, hw={1e3*omega:.3f} meV (FULL Hilbert, no projection)")
    lk0, sp0 = floquet_full(D_op, 0.0, omega)
    print(f"ANCHOR A=0: leak={lk0:.2e}, split={1e3*sp0:.5f} meV vs static {1e3*d0:.5f}")
    print(f"\n{'A(meV)':>7} {'LEAKAGE':>10} {'split(meV)':>11}  note")
    for A in (0.0,0.005,0.010,0.015,0.020,0.030,0.050):
        lk, sp = floquet_full(D_op, A, omega)
        note = "2-level OK" if lk<0.05 else ("MARGINAL" if lk<0.2 else "2-LEVEL BROKEN")
        print(f"{1e3*A:>7.0f} {lk:>10.4f} {1e3*sp:>11.4f}  {note}")
    print("\nnstep convergence (A=30 meV): leakage")
    for nst in (200,300,500):
        lk,_ = floquet_full(D_op,0.030,omega,nstep=nst)
        print(f"   nstep={nst}: leak={lk:.4f}")
