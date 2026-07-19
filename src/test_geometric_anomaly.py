"""Attempt the GEOMETRIC/ANOMALY derivation of a0 proportional to H (the reframe's keystone).
Three linked sub-claims, each with a runnable check:

  A. The BK exponent P(X)~X^{3/2} <=> static phonon energy Integral|grad phi|^3 is SCALE-INVARIANT
     in d=3 spatial dims (the critical/conformal p-Laplacian, p=d). If true, 3/2 is NOT tuned.
  B. Scale invariance FORCES a0 ~ H with a PROTECTED, computable coefficient (the 'anomaly => 1/2pi'
     hope). Test the three candidate scale-breaking carriers of H: non-minimal xi*R; the 4D trace
     anomaly; the Gibbons-Hawking 2pi.
  C. => a0(z) ~ H(z) is a derivation, not an imported hypothesis.

Verdict is decided by the numbers below, not by narrative.
"""
import numpy as np

# ---- units (consistent with the other experiments) ----
hbarc=1.973269804e-5; cm=1.0/hbarc; s=2.99792458e10*cm
Mpc=3.0856775814913673e21*1e3*cm; Mpl=2.435e27
H0=67.0*(1e5*cm)/s/Mpc                       # [eV]
OmDE=0.689; rho_DE=3*Mpl**2*H0**2*OmDE       # [eV^4]
m_len=1.0/1.973269804e-7; s_time=2.99792458e8*m_len; ms2=m_len/s_time**2
a0_obs=1.2e-10*ms2                            # observed a0 [eV]

print("="*76)
print("A. Is P(X)~X^{3/2} <=> Integral |grad phi|^3 scale-invariant in d=3?  (exponent not tuned?)")
print("="*76)
# Build a smooth test field on a 3D grid; compute E_p=Integral|grad phi|^p; rescale phi_lambda(x)=phi(x/L)
# and check E_p[phi_L]/E_p[phi] = L^{d-p}. For p=d=3 the ratio must be 1 (invariant).
N=160; box=16.0; xs=np.linspace(-box/2,box/2,N); h=xs[1]-xs[0]
X,Y,Z=np.meshgrid(xs,xs,xs,indexing='ij')
def base(x,y,z):                              # generic, offset, anisotropic (fixed shape)
    return np.exp(-((x-0.7)**2+y**2+z**2)) + 0.6*np.exp(-(((x+0.9)**2+(y-0.5)**2+z**2)/1.3**2))
def Ep(phi,p):
    gx,gy,gz=np.gradient(phi,h,edge_order=2); g=np.sqrt(gx**2+gy**2+gz**2)
    return np.sum(g**p)*h**3
L=1.5
phi =base(X,   Y,   Z)                         # phi(x)
phiL=base(X/L, Y/L, Z/L)                        # true dilation phi_L(x)=phi(x/L) (centers+widths scaled)
print(f"  {'p':>2} {'E_p[phi]':>12} {'E_p[phi_L]':>12} {'ratio':>8} {'L^(3-p)':>9}  invariant?")
for p in (2,3,4):
    r=Ep(phiL,p)/Ep(phi,p); pred=L**(3-p)
    tag=("YES (scale-free)" if p==3 else f"scales as L^{3-p}")
    ok="  <-- match" if abs(r/pred-1)<0.02 else ""
    print(f"  {p:>2} {Ep(phi,p):12.4f} {Ep(phiL,p):12.4f} {r:8.3f} {pred:9.3f}  {tag}{ok}")
print("  => p=3 (=d) is the UNIQUE scale-invariant power. P(X)~X^{3/2} is the conformal 3-Laplacian,")
print("     singled out by the spatial dimension, NOT fitted.  [A holds]  This is Milgrom (2009),")
print("     'the MOND limit from spacetime scale invariance' (arXiv:0810.4065) -- literature-first.")

print()
print("="*76)
print("B. Does scale invariance FORCE a0 ~ H with a PROTECTED coefficient?  (the load-bearing leap)")
print("="*76)
print(f"  target: a0_obs/H0 = {a0_obs/H0:.4f}   (= 1/2pi = {1/(2*np.pi):.4f}; Milgrom's coincidence)")

print("\n  B1 non-minimal coupling a0^2 = xi*R,  de Sitter R=12 H^2  =>  a0=sqrt(12 xi) H:")
xi_need=(a0_obs/H0)**2/12.0
print(f"     xi needed = (a0/H0)^2/12 = {xi_need:.5f}   vs protected conformal xi=1/6={1/6:.5f}"
      f"  -> off by {(1/6)/xi_need:.0f}x; and xi RUNS (not protected). [B1 fails]")

print("\n  B2 4D trace/Weyl anomaly in de Sitter: <T^mu_mu> ~ (a_anom) H^4  (an ENERGY DENSITY):")
Tanom=H0**4/(2880*np.pi**2)                    # order-of-magnitude de Sitter anomaly density
print(f"     <T> ~ H0^4/(2880 pi^2) = {Tanom:.2e} eV^4   vs rho_DE={rho_DE:.2e} eV^4"
      f"  -> {np.log10(rho_DE/Tanom):.0f} orders too small,")
print(f"     and dimension is [eV^4] (density), NOT [eV] (acceleration). a0 is not an anomaly")
print(f"     coefficient: identifying it with H^4 is the same category error the paper retracted")
print(f"     for 'Lambda^4=condensation free energy'. [B2 fails: category mismatch]")

print("\n  B3 Gibbons-Hawking 2pi: T_dS=H/2pi is a Euclidean-time PERIODICITY factor (thermodynamic),")
print("     not a Weyl-anomaly coefficient. a0=c*T_dS route = Verlinde/Ho-Minic-Ng horizon entropy")
print("     ([11],[13]) -- already known, NOT new, and Test 1 showed the condensate cannot thermally")
print("     couple to that bath (omega/T_dS~1e4). [B3 reduces to Verlinde + is strained by Test 1]")

print()
print("="*76)
print("C. VERDICT")
print("="*76)
print("  A (exponent protected by scale invariance)                          : HOLDS  -> explains")
print("    RAR tightness natively and non-tuned  (== Milgrom 2009).")
print("  B (scale symmetry forces a0 ~ H with computable 1/2pi)              : FAILS  -> the three")
print("    carriers of H all break: xi_R unprotected+wrong number; trace anomaly is a density H^4")
print("    (category error); GH-2pi = Verlinde's horizon route, not an anomaly, and Test-1-strained.")
print("  C (a0(z)~H(z) is a DERIVATION)                                      : NOT ESTABLISHED.")
print("    a0 is the scale where deep-MOND scale-invariance is BROKEN; that scale is set by the")
print("    LOCAL terms in X (DM chemical potential mu(rho), baryon coupling) -- microphysics, not H.")
print("    dS symmetry does not force it to H because mu provides a competing local scale.")
print()
print("  => Keystone does NOT close. Honest residual: a0 proportional to H(z) is an IMPORTED")
print("     hypothesis (Milgrom-1999/Verlinde), testable at high z, NOT a symmetry output. The one")
print("     genuinely protected, non-tuned piece is the EXPONENT (tightness), which is Milgrom 2009.")
