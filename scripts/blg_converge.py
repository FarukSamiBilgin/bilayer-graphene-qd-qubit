"""
Convergence check of the KEY number: |dDelta01/dD_ext| at the d=60nm, eps=0
double-well sweet spot. Vary box size L_box and plane-wave cutoff qmax.
Confinement via ANALYTIC Fourier transform of the double-Gaussian (no Ng grid),
which also cross-checks the earlier FFT-grid result.

Analytic FT: well -U0 exp(-(r-r0)^2/2L^2) -> V(q) = -U0 2pi L^2/Area
             * exp(-q^2 L^2/2) * exp(-i q.r0).
depth mode: r0=(-d/2,0) weight (1+e/2); r0=(+d/2,0) weight (1-e/2).
"""
import numpy as np
from scipy.linalg import eigh
from blg_kp import Hkp, bulk_cond_min

U0, Ldot, d, D_op, h = 0.020, 10.0, 60.0, 0.167, 0.004

def build_basis(L_box, qmax):
    dq=2*np.pi/L_box
    ns=np.arange(-int(qmax/dq)-1,int(qmax/dq)+2)
    QX,QY=[],[]
    for i in ns:
        for j in ns:
            qx,qy=i*dq,j*dq
            if qx*qx+qy*qy<qmax*qmax: QX.append(qx);QY.append(qy)
    return np.array(QX),np.array(QY),L_box

def Vanalytic(qx,qy,L_box,eps):
    pref=-U0*2*np.pi*Ldot**2/L_box**2*np.exp(-(qx**2+qy**2)*Ldot**2/2)
    return pref*((1+eps/2)*np.exp(+1j*qx*d/2)+(1-eps/2)*np.exp(-1j*qx*d/2))

def split(D,eps,QX,QY,L_box):
    n=len(QX); dim=4*n
    H=np.zeros((dim,dim),dtype=complex)
    for a in range(n): H[4*a:4*a+4,4*a:4*a+4]=Hkp(QX[a],QY[a],D,+1)
    for a in range(n):
        for b in range(n):
            val=Vanalytic(QX[a]-QX[b],QY[a]-QY[b],L_box,eps)
            if val!=0:
                for o in range(4): H[4*a+o,4*b+o]+=val
    E=eigh(H,eigvals_only=True)
    cmin=bulk_cond_min(D)
    bnd=E[(E>1e-6)&(E<cmin)]
    return (bnd[1]-bnd[0]) if len(bnd)>=2 else np.nan

def slope_and_split(L_box,qmax):
    QX,QY,Lb=build_basis(L_box,qmax)
    d01=split(D_op,0.0,QX,QY,Lb)
    sl=abs((split(D_op+h,0.0,QX,QY,Lb)-split(D_op-h,0.0,QX,QY,Lb))/(2*h))
    return len(QX),d01,sl

print("Convergence of |dDelta01/dD_ext| at d=60nm sweet spot (analytic confinement)")
print(f"{'L_box':>6}{'qmax':>6}{'N_pw':>6}{'Delta01(meV)':>13}{'|slope|':>10}")
# baseline + L_box variation (qmax=0.5)
for Lb in (100.,120.,140.,160.):
    n,s0,sl=slope_and_split(Lb,0.5)
    print(f"{Lb:>6.0f}{0.5:>6.2f}{n:>6d}{1e3*s0:>13.3f}{sl:>10.4f}")

print("\nBox convergence at SMALLER d (stronger binding). d=30 and d=40nm, qmax=0.5")
print(f"{'d':>4}{'L_box':>7}{'N_pw':>6}{'Delta01(meV)':>13}{'|dDext|':>9}{'d2/deps2':>10}")
heps=0.02
for d_test in (30.0, 40.0):
    globals()['d']=d_test
    for Lb in (120.,150.,180.):
        QX,QY,Lbb=build_basis(Lb,0.5)
        d01=split(D_op,0.0,QX,QY,Lbb)
        sl=abs((split(D_op+h,0.0,QX,QY,Lbb)-split(D_op-h,0.0,QX,QY,Lbb))/(2*h))
        sp=split(D_op,+heps,QX,QY,Lbb); sm=split(D_op,-heps,QX,QY,Lbb)
        d1e=(sp-sm)/(2*heps); d2e=(sp-2*d01+sm)/heps**2
        print(f"{d_test:>4.0f}{Lb:>7.0f}{len(QX):>6d}{1e3*d01:>13.3f}{sl:>9.4f}{d2e:>10.4f}"
              f"  d1e={d1e:.1e}")
