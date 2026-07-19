"""PROPOSAL 1 -- a0 from a COSMIC-UNIFORM density (the tightness-legal reservoir).

Diagnosis: RAR tightness kills a0 that depends on a LOCAL (per-galaxy)
quantity (delta a0/a0 <~ 10-20% allowed; local reservoir densities vary by orders -> dead).
It is BLIND to a0 built from a quantity uniform across galaxies at fixed epoch (delta a0/a0 = 0).

So the only tightness-legal 'reservoir' is a COSMIC density. The natural acceleration from a
density rho is  a ~ c*sqrt(G*rho)  (dimensions: [c][1/time] = accel, since sqrt(G rho)=1/time).
Question: which cosmic density, if any, lands on a0_obs=1.2e-10 m/s^2 -- and what a0(z) does each imply?

This is Door 1 of the safe-zone map. PASS tightness is automatic; the test is (i) does the z=0
AMPLITUDE match without a fudge, and (ii) is the implied O(1) factor DERIVED or just relabeled
(a numerology trap to avoid). Verdict decided by the numbers.
"""
import numpy as np

# ---- constants (SI) ----
c   = 2.99792458e8            # m/s
G   = 6.67430e-11            # m^3 kg^-1 s^-2
Mpc = 3.0856775814913673e22  # m
H0  = 67.0e3/Mpc             # s^-1
a0  = 1.2e-10               # observed MOND scale [m/s^2]

Om, OL, Ob = 0.31, 0.69, 0.049
rho_crit = 3*H0**2/(8*np.pi*G)          # kg/m^3
DENS = [
    ("baryons        rho_b = Om_b rho_c", Ob*rho_crit, 3.0),   # (name, rho, d ln rho/d ln(1+z))
    ("matter         rho_m = Om_m rho_c", Om*rho_crit, 3.0),
    ("dark energy    rho_L = Om_L rho_c", OL*rho_crit, 0.0),
    ("critical       rho_crit (~ H^2)",   rho_crit,    None),   # ~H^2 -> a ~ H (rate branch)
]

print("="*82)
print("PROPOSAL 1: a0 = c*sqrt(G*rho) for cosmic-uniform densities  (tightness-legal by design)")
print("="*82)
print(f"  a0_obs = {a0:.2e} m/s^2 ;  c H0 = {c*H0:.2e} ;  a0/(cH0) = {a0/(c*H0):.4f}")
print(f"  Milgrom/Verlinde coincidence 1/2pi = {1/(2*np.pi):.4f}  (=> cH0/2pi is {c*H0/(2*np.pi)/a0:.2f}*a0,"
      f" the famous 16% low)")
print()
print(f"  {'cosmic density':<34}{'c*sqrt(G rho) [m/s^2]':>22}{'/ a0_obs':>10}{'implied O(1)=a/cH0':>20}")
for name, rho, _ in DENS:
    a = c*np.sqrt(G*rho)
    print(f"  {name:<34}{a:>22.3e}{a/a0:>10.3f}{a/(c*H0):>20.4f}")

print()
print("="*82)
print("The implied a0(z) branch for each  (a ~ sqrt(rho); rho(z) known)")
print("="*82)
def Ez(z): return np.sqrt(Om*(1+z)**3 + OL)
zs = [0.5, 1, 2, 3]
print(f"  {'branch':<30}" + "".join(f"z={z:<7}" for z in zs))
# sqrt(rho_L)=const; sqrt(rho_crit)~H -> E(z); sqrt(rho_m or rho_b)~(1+z)^1.5
rows = [
    ("a0 ~ sqrt(rho_L)  (Milgrom)",  lambda z: 1.0+0*z),
    ("a0 ~ sqrt(rho_crit) ~ H (Verlinde)", lambda z: Ez(z)),
    ("a0 ~ sqrt(rho_m)  (=(1+z)^1.5)",     lambda z: (1+z)**1.5),
]
for name, f in rows:
    print(f"  {name:<30}" + "".join(f"{f(z):<9.2f}" for z in zs))

print()
print("="*82)
print("VERDICT")
print("="*82)
a_m = c*np.sqrt(G*Om*rho_crit)
print(f"  (i) TIGHTNESS: PASS for all -- each rho is uniform across galaxies at fixed z, so")
print(f"      delta a0/a0 = 0 galaxy-to-galaxy. This is what the LOCAL reservoir branch could not do.")
print(f"  (ii) AMPLITUDE: matter density gives c*sqrt(G rho_m) = {a_m:.2e} = {a_m/a0:.2f}*a0 -- a 4% match,")
print(f"      BETTER than the 2pi coincidence (16%). rho_L -> 1.56*a0, rho_crit -> 1.87*a0, rho_b -> 0.42*a0.")
print(f"  (iii) DERIVED? NO. a/cH0 = sqrt(3 Om/8pi) = {np.sqrt(3*Om/(8*np.pi)):.3f} is an undetermined O(1),")
print(f"      exactly like 1/2pi = {1/(2*np.pi):.3f}. Both are 'a0 ~ cH0 * O(1)'. The 4% is a better-fitting")
print(f"      free number, NOT a derivation -> this is a numerology trap, tightness-legal")
print(f"      but not explanatory. Proposal 1 JOINS the graveyard AS A MECHANISM.")
print()
print(f"  USEFUL BY-PRODUCT (kept): writing a0 ~ sqrt(rho_X) makes each choice a DISTINCT, tightness-legal")
print(f"  a0(z). The z=0 amplitude marginally prefers rho_MATTER (1.04) over rho_Lambda (1.56) -- a")
print(f"  physically different story (a0 from gravitating matter, Mach/Sciama-like, not from vacuum) and")
print(f"  the STEEPEST evolution: a0(z=2) ~ (1+2)^1.5 = {(3)**1.5:.1f}x, vs rate branch E(2)={Ez(2):.1f}x,")
print(f"  vs Lambda ~1x. This adds a 4th curve to the a0(z) shape test and is the one live output here.")
print(f"  CAVEAT: a0~c*sqrt(G rho_m) is very likely a rearrangement of known coincidences -- LIT-CHECK")
print(f"  before any novelty claim (the repo's own #1 lesson).")
