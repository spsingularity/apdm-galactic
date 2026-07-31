#!/usr/bin/env python3
"""
E-1' (proxy) : Does the rising-a0 trend survive the mass confound, per-galaxy?

E-1 used the PUBLISHED MUSE global fit a0(z) = 1.00 + 1.59z, which carries the
astrophysical transfer B(z) (mass/selection evolution) inside it. The released
per-galaxy catalogue (muse_release_galaxies.csv: z, log10 M*, median outer rotational
acceleration a_rot for r > 2 kpc) lets us CONTROL the dominant confound directly:

    log10 a_rot = alpha + beta * z + m * log10 M*      (multiple regression)

In the deep-MOND regime a_rot ~ sqrt(a0 g_bar), so at fixed mass (g_bar fixed up to
size evolution) beta = 0.5 * dlog10 a0/dz. Predictions over the MUSE window:
    rate branch (a0 ~ H)      : beta ~ +0.11
    flat                      : beta =  0
    density (a0 ~ sqrt(rhoDE)): beta ~ -0.02
    MUSE global fit           : beta ~ +0.14
This is a CONSISTENCY PROBE, not a measurement of a0 (no per-galaxy g_bar in the
release; radius coverage and surface-brightness selection can still bias beta).
"""
import numpy as np, csv, os

HERE = os.path.dirname(__file__)
rows = list(csv.DictReader(open(os.path.join(HERE, "..", "results", "muse_release_galaxies.csv"))))

# numeric candidates with the paper's numeric cuts (mass, z, logZ), a_rot available
sel = [r for r in rows if r["numeric_candidate"]=="True" and float(r["median_arot_r_gt_2kpc_m_s2"])>0]
z    = np.array([float(r["z"]) for r in sel])
lm   = np.array([float(r["log10_mstar_msun"]) for r in sel])
lar  = np.log10([float(r["median_arot_r_gt_2kpc_m_s2"]) for r in sel])

print("="*78)
print("E-1' proxy: mass-controlled per-galaxy a_rot(z) trend (MUSE release)")
print("="*78)
print(f"  N = {len(sel)} numeric candidates (paper's final sample: 79 of these 85)")
print(f"  z range [{z.min():.2f}, {z.max():.2f}]   log10 M* range [{lm.min():.2f}, {lm.max():.2f}]")

# the confound itself: does the sample get more massive with z?
cz = np.corrcoef(z, lm)[0,1]
print(f"  mass-redshift selection correlation corr(log M*, z) = {cz:+.2f}"
      f"   -> the confound is {'REAL and must be controlled' if abs(cz)>0.15 else 'mild'}")

def ols(X, y):
    XtX = X.T@X; b = np.linalg.solve(XtX, X.T@y)
    res = y - X@b; s2 = res@res/(len(y)-X.shape[1])
    return b, np.sqrt(np.diag(s2*np.linalg.inv(XtX)))

# (a) raw trend (no mass control)
b, e = ols(np.column_stack([np.ones_like(z), z]), lar)
print(f"\n  (a) raw:            dlog10 a_rot/dz = {b[1]:+.3f} +/- {e[1]:.3f}")
# (b) mass-controlled
X = np.column_stack([np.ones_like(z), z, lm - lm.mean()])
b2, e2 = ols(X, lar)
print(f"  (b) mass-controlled: dlog10 a_rot/dz = {b2[1]:+.3f} +/- {e2[1]:.3f}"
      f"   (mass slope m = {b2[2]:+.2f} +/- {e2[2]:.2f})")

print(f"\n  branch predictions for beta = 0.5 dlog10 a0/dz over this window:")
for lab, beta in [("rate  (a0 ~ H)", 0.5*np.log10(1.888)/1.11),
                  ("flat", 0.0),
                  ("density (a0 ~ sqrt rho, CPL)", 0.5*np.log10(0.894)/1.11),
                  ("MUSE global fit (1+1.59z)", 0.5*1.59/(1+1.59*z.mean())/np.log(10))]:
    t = (b2[1]-beta)/e2[1]
    print(f"    {lab:30}: beta = {beta:+.3f}   -> {t:+.1f} sigma from measured (b)")

print(f"""
  READING (state honestly):
   * A positive mass-controlled beta consistent with ~+0.11 supports a genuinely
     rising a0(z) (rate branch) beyond the mass confound.
   * beta consistent with 0 after mass control would mean the published rising a0
     is largely selection/transfer B(z) -- undercutting E-1's rate-branch support.
   * a_rot is NOT the RAR-fitted a0: per-galaxy g_bar, size evolution, and pressure
     support are uncontrolled here. Treat as a directional consistency check; the
     full E-1' requires the unreleased baryonic decompositions.
""")
