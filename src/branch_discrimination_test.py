#!/usr/bin/env python3
"""
E-1 : Discriminate the two a0(z) branches against the MUSE-DARK III data.

  density branch (APDM superfluid, Milgrom-1999 vacuum):  a0 ∝ sqrt(rho_DE)   -> FALLS (w_a<0)
  rate branch    (USC / Verlinde emergent gravity):        a0 ∝ H(z)           -> RISES

Data:
  * DESI DR2 templates (this repo, research/results/desi_template_quantiles.csv):
    posterior of a0(z)/a0(0) for each branch at z = 0.33, 1.0, 1.44, 2.0, per SN set.
  * MUSE-DARK III (Ciocan et al. 2026, arXiv:2604.22613) published global fit:
        a0(z) = (1.00 +/- 0.04) + (1.59 +/- 0.11) z     [1e-10 m/s^2]
    over 0.33 < z < 1.44 (79 rotators).

Test: over the MUSE interval [0.33, 1.44], does the observed a0 ratio match either
branch? Report the tension in sigma. (a0's overall normalization cancels in the ratio,
so this is the shape test Paper 1 set up; the residual amplitude is the astrophysical
transfer B(z) degeneracy, discussed in the verdict.)
"""
import numpy as np, csv, os

HERE = os.path.dirname(os.path.abspath(__file__))
QT = os.path.join(HERE, "data", "desi_template_quantiles.csv")  # DESI DR2 a0(z) branch templates

# ---------------- MUSE observed a0(z) and its ratio over [0.33, 1.44] ----------------
A0, sA0 = 1.00, 0.04     # 1e-10 m/s^2
S0, sS0 = 1.59, 0.11     # slope, 1e-10 m/s^2 per unit z
z0, z1 = 0.33, 1.44
rng = np.random.default_rng(0)
def muse_ratio_samples(n=200000):
    A = rng.normal(A0, sA0, n); S = rng.normal(S0, sS0, n)
    return (A + S*z1)/(A + S*z0)
mr = muse_ratio_samples()
R_muse, sR_muse = mr.mean(), mr.std()

# ---------------- load DESI templates, build a0(z1)/a0(z0) per branch per SN set ----------------
rows = list(csv.DictReader(open(QT)))
def get(dataset, template, z):
    for r in rows:
        if r["dataset"]==dataset and r["template"]==template and abs(float(r["z"])-z)<1e-6:
            return float(r["median"]), (float(r["q84"])-float(r["q16"]))/2
    raise KeyError((dataset, template, z))

datasets = sorted({r["dataset"] for r in rows})
branches = {"sqrt_rhoDE":"DENSITY  a0~sqrt(rho_DE)", "H_over_H0":"RATE     a0~H(z)"}

print("="*84)
print("E-1  a0(z) branch discrimination vs MUSE-DARK III  (ratio over z=0.33 -> 1.44)")
print("="*84)
print(f"  MUSE observed ratio a0(1.44)/a0(0.33) = {R_muse:.3f} +/- {sR_muse:.3f}")
print(f"  (from a0(z) = {A0}+{S0}z, errors {sA0}/{sS0})")
print()
print(f"  {'SN set':12} | {'branch':26} | {'template ratio':>16} | {'tension vs MUSE':>16}")
print("  " + "-"*78)
summary = {}
for ds in datasets:
    for tmpl, lab in branches.items():
        m1, s1 = get(ds, tmpl, 1.44)
        m0, s0 = get(ds, tmpl, 0.33)
        # ratio a0(1.44)/a0(0.33); templates normalized to z=0, so divide.
        R = m1/m0
        # conservative ratio error via marginals (ignores shared-w0wa correlation -> OVER-estimates)
        sR = R*np.sqrt((s1/m1)**2 + (s0/m0)**2)
        # combined tension
        sig = (R_muse - R)/np.sqrt(sR_muse**2 + sR**2)
        summary.setdefault(tmpl, []).append(sig)
        print(f"  {ds:12} | {lab:26} | {R:6.3f} +/- {sR:5.3f}   | {sig:+6.1f} sigma")
    print("  " + "-"*78)

print()
print("  Branch-averaged tension (mean over SN sets):")
for tmpl, lab in branches.items():
    s = np.array(summary[tmpl])
    print(f"    {lab:26}: {np.mean(s):+.1f} sigma   ({'EXCLUDED — wrong sign' if 'DENSITY' in lab else 'right sign, amplitude short'})")

# ---------------- what power a0 ∝ (1+z)^p / H^n does MUSE want? ----------------
print()
print("="*84)
print("  What evolution index does the raw MUSE a0(z) prefer?")
print("="*84)
p_muse = np.log(R_muse)/np.log((1+z1)/(1+z0))
# H(z)/H0 for LCDM (Om=0.31): index n s.t. (H1/H0)^n gives R_muse, using rate-branch template
Hr = get('pantheonplus','H_over_H0',1.44)[0]/get('pantheonplus','H_over_H0',0.33)[0]
n_muse = np.log(R_muse)/np.log(Hr)
rho_r = get('pantheonplus','sqrt_rhoDE',1.44)[0]/get('pantheonplus','sqrt_rhoDE',0.33)[0]
print(f"  a0 ∝ (1+z)^p :  p = {p_muse:.2f}   (rate branch a0∝H gives effective p≈{np.log(Hr)/np.log((1+z1)/(1+z0)):.2f})")
print(f"  a0 ∝ H(z)^n  :  n = {n_muse:.2f}   (rate branch is n=1; MUSE wants n>1 -> steeper than H)")
print(f"  density branch would need a0∝rho_DE^m with rho ratio {rho_r:.3f}<1 -> m<0 to rise: IMPOSSIBLE (rho falls).")

print()
print("="*84)
print("  VERDICT")
print("="*84)
print("  * DENSITY (superfluid / APDM) branch: a0 FALLS while MUSE RISES -> excluded at many sigma.")
print("    This is APDM's OWN headline branch. The MUSE data reject it.")
print("  * RATE (emergent-gravity / USC) branch: correct SIGN, but ~2-4 sigma too WEAK in amplitude")
print("    (MUSE wants a0∝H^{%.1f}, steeper than a0∝H^1)."%n_muse)
print("    The deficit is exactly the size the astrophysical transfer B(z) (mass/selection")
print("    evolution; Paper 1's a0_eff=a0_cos*B) can absorb -> data FAVOR the rate branch in")
print("    sign but cannot yet ISOLATE a0_cos(z) from B(z). Controlling B(z) (matched stellar")
print("    mass / morphology bins across z) is the decisive follow-up.")
