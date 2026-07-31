#!/usr/bin/env python3
"""D-4' execution: the accelerated-worldline influence-functional kernel for the
T^{mu nu} grad_mu X_nu channel, in the conformal-scalar proxy model.

Structure derived analytically (verified numerically here):
 (1) On a worldline, int T grad X = -m int dtau a^nu X_nu  (since div T = m a^nu delta^4):
     the displacement field couples DIRECTLY to the acceleration vector.
 (2) Second-order influence functional on the uniformly accelerated (Deser-Levin)
     trajectory in the dS Bunch-Davies state: the pulled-back Wightman function is
     thermal with kappa = sqrt(a^2 + H^2)  (Deser-Levin theorem):
        G+(dtau) = -(kappa^2/16 pi^2) / sinh^2(kappa(dtau - i eps)/2)
     and the acceleration two-point contraction on the flat-embedding trajectory is
        a^nu(tau) a_nu(tau') = a^2 cosh(a dtau).
 (3) The kernel integral (exact):
        K(a,H) = int dtau cosh(a dtau) G+(dtau) = (a / 4pi) * cot(pi a / kappa)
     via  int dx cosh(alpha x)/sinh^2(beta(x-i eps)) = -(pi alpha/beta^2) cot(pi alpha/(2 beta)),
     convergent iff alpha < 2 beta  <=>  a < kappa  (always true; marginal as H->0).

This script verifies the master integral numerically and maps the physics of K(a,H).
"""
import numpy as np
from scipy.integrate import quad

# ---------- 1. verify the master integral ----------
def master_numeric(alpha, beta, L=60.0):
    # Deform the contour to Im x = -pi/(2 beta) (crosses no poles; the strip between the
    # pole rows at Im=0 and Im=-pi/beta is analytic). There sinh(beta(x)) -> -i cosh(beta t),
    # so the integrand is smooth: cosh(alpha(t - i pi/(2beta)))/(-cosh^2(beta t)).
    # The odd (imaginary) part integrates to zero; the real part is
    #   -cos(pi alpha/(2 beta)) * int dt cosh(alpha t)/cosh^2(beta t).
    def f(t):
        return np.cosh(alpha*t)/np.cosh(beta*t)**2
    val,_ = quad(f, -L, L, limit=400)
    return -np.cos(np.pi*alpha/(2*beta))*val

def master_closed(alpha, beta):
    return -(np.pi*alpha/beta**2)/np.tan(np.pi*alpha/(2*beta))

print("="*74)
print("D-4' worldline kernel: verify  int cosh(ax)/sinh^2(b(x-ie)) = -(pi a/b^2) cot(pi a/2b)")
print("="*74)
print("  (numeric via pole-free contour Im x = -pi/2beta; equals -cos(pi a/2b) * int cosh(at)/cosh^2(bt),")
print("   which by Gradshteyn 3.512 = -(pi a/b^2) cot(pi a/2b) -- i.e. the closed form is ALSO proven analytically)")
for alpha,beta in [(0.3,1.0),(0.9,1.0),(1.5,1.0),(0.5,0.7),(1.2,0.8)]:
    num=master_numeric(alpha,beta); clo=master_closed(alpha,beta)
    print(f"  alpha={alpha}, beta={beta}:  numeric={num:+.6f}   closed={clo:+.6f}   rel.err={abs(num-clo)/abs(clo):.1e}")

# ---------- 2. the physical kernel K(a,H) ----------
# K(a,H) = (m~^2 a^2/2) * int cosh(a dtau) G+(dtau)  with  G+ = -(kappa^2/16pi^2) sinh^-2(kappa dtau/2)
#        = (m~^2 a^2/2) * (kappa^2/16pi^2) * (pi a/(kappa/2)^2) cot(pi a/kappa)   [note sign: -G+ * -(...)]
#        = (m~^2/2) * (a^3/(4 pi)) cot(pi a / kappa)
def K(a, H=1.0):
    kap = np.sqrt(a*a + H*H)
    return (a**3/(4*np.pi))/np.tan(np.pi*a/kap)

print()
print("="*74)
print("K(a,H)/m~^2 = (a^3/8pi) cot(pi a/kappa),  kappa = sqrt(a^2+H^2)   [H=1 units]")
print("="*74)
print(f"  {'a/H':>8} | {'K':>12} | {'K_deepMOND = a^2 H/(8 pi^2)':>26} | note")
for a in [0.01,0.05,0.1,0.3,1/np.sqrt(3),0.6,1.0,2.0,5.0]:
    deep = a*a/(8*np.pi**2)   # leading small-a limit: (a^3/8pi)*(kappa/pi a) -> a^2 kappa/(8 pi^2) -> a^2 H/8pi^2
    note = "<-- ZERO of kernel at a = H/sqrt(3)" if abs(a-1/np.sqrt(3))<1e-9 else ""
    print(f"  {a:8.3f} | {K(a)/2:+12.5e} | {deep/1:26.5e} | {note}")

print()
print("limits:")
print("  a << H : K -> m~^2 a^2 H / (16 pi^2)      [deep regime: quadratic in a, linear in H]")
print("  a -> H/sqrt(3): K = 0                     [kernel zero: pi a/kappa = pi/3... check]")
# check where cot(pi a/kappa)=0: pi a/kappa = pi/2 => a = kappa/2 => 4a^2 = a^2+H^2 => a = H/sqrt(3)
a0 = 1/np.sqrt(3)
print(f"  cot zero: pi*a/kappa = pi/2 at a/H = 1/sqrt(3) = {a0:.4f}: K({a0:.4f}) = {K(a0):.2e}")
print("  a >> H : K -> -(m~^2/2) a^5/(4 pi^2 H^2)  [DIVERGES as H->0: eternal-Rindler IR pathology")
print("           regularized only by H; signals the idealization breaks, finite worldline needed]")

# ---------- 3. compare crossover scales ----------
print()
print("="*74)
print("crossover-scale comparison (H0 = 70 km/s/Mpc, c=1 units internal)")
print("="*74)
cH0 = 6.80e-10  # m/s^2 for H0=70
a0_obs = 1.2e-10
print(f"  observed a0            = {a0_obs:.2e} m/s^2 = {a0_obs/cH0:.4f} cH0")
print(f"  conjecture cH0/2pi     = {cH0/(2*np.pi):.2e} m/s^2 = {1/(2*np.pi):.4f} cH0   (obs/conj = {a0_obs/(cH0/(2*np.pi)):.3f})")
print(f"  kernel zero H/sqrt(3)  = {cH0/np.sqrt(3):.2e} m/s^2 = {1/np.sqrt(3):.4f} cH0   (obs/zero = {a0_obs/(cH0/np.sqrt(3)):.3f})")
print(f"  raw Deser-Levin 2cH    = {2*cH0:.2e} m/s^2 = 2.0 cH0                (obs/DL   = {a0_obs/(2*cH0):.3f})")
print()
print("  NOTE the kernel normalization carries EXACTLY one 1/(4pi):")
print("     int cosh(a t) G+(t) dt = (a/4pi) cot(pi a/kappa)")
print("  -> the conjectured Green-function 4pi DOES appear at the conjectured place;")
print("     but the interpolation is cot(pi a/kappa), NOT Milgrom's sqrt(a^2+H^2)-H,")
print("     and the kernel VANISHES at a=H/sqrt(3) (so it cannot literally BE the inertia).")
