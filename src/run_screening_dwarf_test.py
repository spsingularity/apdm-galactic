"""Field-vs-satellite-dwarf SCREENING test (falsifier for a candidate mechanism, the Phi-screened
DM-DE scalar).

THE DISCRIMINANT. A chameleon/symmetron-type scalar is screened by POTENTIAL DEPTH Phi.
Plain MOND's environmental suppression (the External Field Effect, EFE) is set by the
EXTERNAL ACCELERATION g_ext = |grad Phi|. These are DIFFERENT variables, and they
DECORRELATE: a dwarf near a cluster CENTRE has a deep potential Phi (~ v_cluster^2/c^2)
but a SMALL local acceleration g_ext (-> 0 at the centre). There the two hypotheses make
OPPOSITE predictions for the radial-acceleration-relation (RAR) residual:

  Phi-screening : deep Phi -> SCREENED -> extra force off -> Newtonian (RAR residual -> 0)
  MOND-EFE      : small g_ext -> UNSCREENED -> full MOND boost (RAR residual -> full)
  no-EFE mod-gravity : full MOND boost regardless of environment
  CDM           : mass discrepancy set by the halo; environment enters only via tidal
                  stripping (tracks tidal field ~ g_ext / pericentre, NOT Phi)

So the RAR residual's PARTIAL correlation with Phi (at fixed g_ext) is the screening
signature; its partial correlation with g_ext (at fixed Phi) is the MOND-EFE signature.
Cluster-central dwarfs are the golden decorrelated sample.

This script is a FORWARD MODEL + DISCRIMINATION FORECAST (it does not fabricate data). To
confront real data: SPARC + Local-Volume dwarf kinematics x a screening map of the local
universe (Phi from 2M++/BORG reconstruction, Desmond & Ferreira 2019). No DM detection needed.

Run: python3 experiments/run_screening_dwarf_test.py
"""
import numpy as np
rng = np.random.default_rng(1)

# ---- scales ----
c_kms = 2.99792458e5
a0 = 1.0                                   # work in units of a0 (=1.2e-10 m/s^2)
Phi_crit = 1e-6                            # screening threshold (Desmond-Ferreira chi_c ~ 1e-7..1e-5)

def Phi(v_kms):                            # potential depth ~ v^2/c^2 (dimensionless)
    return (v_kms/c_kms)**2

# ---- environment classes: (name, Phi_env, g_ext/a0, N_available) ----
# Phi from host/LSS velocity scale; g_ext = local external acceleration (small at cluster centre).
classes = [
    # name                         v_env[km/s]  g_ext/a0   comment
    ("field / isolated dwarf",         30,        0.02,   "shallow Phi, tiny g_ext"),
    ("MW-like satellite",             220,        0.5,    "moderate Phi, moderate g_ext"),
    ("cluster dwarf, mid-radius",    1000,        2.0,    "deep Phi AND large g_ext"),
    ("cluster dwarf, CENTRAL",       1000,        0.03,   "deep Phi, SMALL g_ext  <-- golden"),
]

g_bar = 0.1                                # internal baryonic acceleration (deep-MOND regime), /a0
Delta_full = 0.5*np.log10(a0/g_bar)        # deep-MOND boost = 1/2 log10(a0/g_bar)

def pred(hyp, Phi_env, g_ext):
    """Predicted RAR residual Delta = log10(g_obs/g_bar) for each hypothesis."""
    if hyp == "Phi-screen":
        return 0.0 if Phi_env > Phi_crit else Delta_full          # screened -> Newtonian
    if hyp == "MOND-EFE":
        return Delta_full / (1.0 + g_ext/a0)                       # boost dies as g_ext>a0
    if hyp == "modgrav-noEFE":
        return Delta_full                                          # no environment dependence
    if hyp == "CDM":
        return Delta_full                                          # RAR is emergent; ~env-independent
    raise ValueError(hyp)

print(f"### Predicted RAR residual Delta (dex)   [g_bar={g_bar}a0, Delta_full={Delta_full:.2f}] ###")
hyps = ["Phi-screen", "MOND-EFE", "modgrav-noEFE", "CDM"]
print(f"  {'environment':28s} {'Phi_env':>9s} {'g_ext/a0':>8s} | " + " ".join(f"{h:>13s}" for h in hyps))
for name, v, ge, _ in classes:
    Pe = Phi(v)
    row = " ".join(f"{pred(h,Pe,ge):13.2f}" for h in hyps)
    print(f"  {name:28s} {Pe:9.1e} {ge:8.2f} | {row}")
print("  -> GOLDEN ROW (cluster CENTRAL): Phi-screen says 0.00 (Newtonian), MOND-EFE says"
      f" {pred('MOND-EFE',Phi(1000),0.03):.2f} (near-full MOND). OPPOSITE predictions.\n")

# ---- discrimination power on the golden (cluster-central) sample ----
print("### Discrimination: Phi-screen vs MOND-EFE on cluster-CENTRAL dwarfs ###")
d_screen = pred("Phi-screen", Phi(1000), 0.03)
d_efe    = pred("MOND-EFE",   Phi(1000), 0.03)
gap = abs(d_efe - d_screen)
sig_dex = 0.10                              # per-object RAR-residual scatter (intrinsic ~0.1 dex)
print(f"  |Delta_EFE - Delta_screen| = {gap:.2f} dex ;  per-object scatter = {sig_dex:.2f} dex")
for N in (5, 10, 25, 50):
    signif = gap/(sig_dex/np.sqrt(N))
    print(f"    N_central={N:3d}  ->  separation = {signif:4.1f} sigma")
print("  STATISTICALLY, ~10-25 cluster-central dwarfs with clean internal kinematics decide it.\n")

# ---- partial-correlation structure (the general, non-golden test) ----
print("### Partial-correlation signature (full sample, controlling confounds) ###")
# Build a synthetic population spanning the classes; recover which variable drives Delta.
def population(truth, n=4000):
    picks = rng.integers(0, len(classes), n)
    Pe = np.array([Phi(classes[i][1]) for i in picks]) * rng.lognormal(0, 0.3, n)
    ge = np.array([classes[i][2] for i in picks]) * rng.lognormal(0, 0.3, n)
    D  = np.array([pred(truth, Pe[k], ge[k]) for k in range(n)]) + rng.normal(0, sig_dex, n)
    return np.log10(Pe), np.log10(ge), D
def partial_corr(x, y, z):                  # corr(x,y) controlling for z
    def resid(a, b):
        b1 = np.c_[np.ones_like(b), b]
        return a - b1 @ np.linalg.lstsq(b1, a, rcond=None)[0]
    rx, ry = resid(x, z), resid(y, z)
    return np.corrcoef(rx, ry)[0, 1]
print(f"  {'if truth =':16s} {'pcorr(D,Phi|g_ext)':>18s} {'pcorr(D,g_ext|Phi)':>19s}   -> reads as")
for truth in ("Phi-screen", "MOND-EFE", "modgrav-noEFE"):
    lP, lg, D = population(truth)
    pP = partial_corr(D, lP, lg); pG = partial_corr(D, lg, lP)
    reads = "SCREENING" if abs(pP) > abs(pG)+0.1 else ("MOND-EFE" if abs(pG) > abs(pP)+0.1 else "no env-dep")
    print(f"  {truth:16s} {pP:+18.2f} {pG:+19.2f}   -> {reads}")
print("  -> the residual correlates with Phi (screening) OR g_ext (EFE) -- the partial")
print("     correlations cleanly separate the hypotheses even without the golden sample.\n")

print("### HONEST CAVEATS ###")
print("  * Statistical power is high (~tens of dwarfs); the REAL barrier is SYSTEMATIC:")
print("    clean internal kinematics of cluster-central dwarfs are hard (small, distant, crowded).")
print("  * CDM TIDAL STRIPPING is the main confound -- it also makes satellites deviate, but tracks")
print("    tidal field/pericentre (~g_ext), not Phi; control by excluding morphologically disturbed")
print("    dwarfs and by the partial-correlation-with-Phi-at-fixed-g_ext test.")
print("  * Result depends on Phi_crit (scanned ~1e-7..1e-5); if Phi_crit > cluster Phi (~1e-5),")
print("    NOTHING screens and the test is void -> Phi_crit is itself a prediction to pin/vary.")
print("  * Baryonic RAR scatter (inclination, distance, gas) sets the 0.1 dex floor used above.")
