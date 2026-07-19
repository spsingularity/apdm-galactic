"""PAPER 3A PROTOTYPE -- background-vs-growth consistency null test of the DESI w0wa "hint".

MOTIVATION. At the background level a coupled dark-matter--dark-energy sector is EXACTLY degenerate
with uncoupled evolving dark energy (Kunz, astro-ph/0702615): geometry measures only rho_dm+rho_de.
The degeneracy breaks at perturbation level: real energy/momentum exchange alters DM growth
(G_eff != G, extra drag). So: fit the BACKGROUND to DESI DR2 BAO alone, predict fsigma8(z) under
GR + uncoupled smooth DE with ONE free amplitude sigma8_0, and form the null ratio
    O(z_i) = fsigma8_obs(z_i) / fsigma8_pred(z_i).
Because sigma8_0 is fit to the same points, the inverse-variance weighted mean of O is 1 BY
CONSTRUCTION -- the test lives entirely in the SHAPE of O(z): chi2 of O=1 and any coherent
deviation localized at the redshifts of the DESI w0wa hint (z ~ 0.3-0.8).

PRE-REGISTERED READ (written before any numbers were computed):
  * O(z)=1 within ~2sigma everywhere (chi2/dof of O=1 not anomalous, |band mean - 1| < 2sigma for
    the z=0.3-0.8 band, under BOTH backgrounds) => NO perturbation-level support for a coupled
    sector; the DESI evolving-DE hint stays background-only, where coupling and uncoupled w(z) are
    indistinguishable. Log as a (useful) null.
  * A coherent >2sigma deviation of the z=0.3-0.8 band mean, present under the w0waCDM background
    (i.e. NOT absorbed by the extra background freedom) and consistent between the DESI-DR1-FS and
    classic-eBOSS subsets => genuinely interesting; warrants the full paper with proper covariances.
  * Single-point >2sigma excursions with no band coherence => noise; still a null.
  * Sanity gates (must pass or the verdict is void): w0waCDM fit to DR2 BAO prefers w0 > -1,
    wa < 0; LCDM Om in ~[0.28, 0.32].

DATA (all hardcoded from the source papers; no downloads needed):
  * DESI DR2 BAO distances, arXiv:2503.14738 Table (results section, main.tex l.565-585 of the
    arXiv source): D_V/r_d for BGS; (D_M/r_d, D_H/r_d, r_MH) for LRG1, LRG2, LRG3+ELG1, ELG2,
    QSO, Lya.  [REAL]
  * DESI DR1 full-shape ShapeFit+BAO growth, arXiv:2411.12021 Table 6 (tab:SF-results,
    fsigma_s8/fid ratios, MAP +/- 68%) x Table (tab:fiducial_values, appendix) fiducial
    fsigma_s8(z); fsigma8 = ratio x fiducial (assumes fiducial r_d = true r_d, as the paper's own
    Fig. 4 does; errors symmetrized).  [REAL]
  * Classic compilation: 6dFGS (arXiv:1204.4725); SDSS MGS + BOSS DR12 + eBOSS DR16 consensus
    BAO+RSD row of Alam+2021 arXiv:2007.08991 Table 3.  [REAL]

KNOWN CAVEAT (stated up front, cannot be fixed in a prototype): every fsigma8 point was measured
against a fiducial-cosmology template (Alcock-Paczynski + dilation + fiducial r_d for the DESI
fsigma_s8 -> fsigma8 conversion). A full analysis must propagate the template dependence and the
cross-bin covariances; a prototype treats points as independent and template-corrected.

Run: python3 experiments/test_desi_growth_null.py   (numpy + scipy only)
"""
import numpy as np
from scipy import integrate, optimize, stats

C_KMS = 299792.458

# ----------------------------------------------------------------------------------------------
# 1. DESI DR2 BAO data (arXiv:2503.14738, results table; exact printed values)
#    columns: z_eff, DM/rd, sig, DH/rd, sig, r_MH   (BGS is DV/rd only)
# ----------------------------------------------------------------------------------------------
BAO_DV = [("BGS", 0.295, 7.942, 0.075)]                                   # D_V/r_d
BAO_MH = [  # tracer, z_eff, DM/rd, err, DH/rd, err, corr(DM,DH)
    ("LRG1",      0.510, 13.588, 0.167, 21.863, 0.425, -0.459),
    ("LRG2",      0.706, 17.351, 0.177, 19.455, 0.330, -0.404),
    ("LRG3+ELG1", 0.934, 21.576, 0.152, 17.641, 0.193, -0.416),
    ("ELG2",      1.321, 27.601, 0.318, 14.176, 0.221, -0.434),
    ("QSO",       1.484, 30.512, 0.760, 12.817, 0.516, -0.500),
    ("Lya",       2.330, 38.988, 0.531,  8.632, 0.101, -0.431),
]

# ----------------------------------------------------------------------------------------------
# 2. fsigma8 compilation
# ----------------------------------------------------------------------------------------------
# DESI DR1 full-shape (arXiv:2411.12021): ShapeFit+BAO MAP ratios fsig_s8/fid (Table 6, errors
# symmetrized as mean of +/-) times the fiducial fsig_s8(z) (appendix tab:fiducial_values).
DESI_FS = [  # tracer, z_eff, ratio, sym err, fiducial fsigma_s8
    ("DESI BGS",  0.295, 0.84,  0.19,             0.4723),
    ("DESI LRG1", 0.510, 1.16,  0.13,             0.4733),
    ("DESI LRG2", 0.706, 1.04,  (0.11+0.092)/2,   0.4608),
    ("DESI LRG3", 0.919, 0.997, (0.10+0.084)/2,   0.4398),
    ("DESI ELG2", 1.317, 0.945, (0.097+0.077)/2,  0.3944),
    ("DESI QSO",  1.491, 1.16,  0.12,             0.3750),
]
# Classic: 6dFGS (arXiv:1204.4725); MGS/BOSS/eBOSS = BAO+RSD consensus row, Alam+2021
# arXiv:2007.08991 Table 3 (z=0.15, 0.38, 0.51, 0.70, 0.85, 1.48).
CLASSIC = [
    ("6dFGS",     0.067, 0.423, 0.055),
    ("SDSS MGS",  0.150, 0.530, 0.160),
    ("BOSS z1",   0.380, 0.497, 0.045),
    ("BOSS z2",   0.510, 0.459, 0.038),
    ("eBOSS LRG", 0.700, 0.473, 0.041),
    ("eBOSS ELG", 0.850, 0.315, 0.095),
    ("eBOSS QSO", 1.480, 0.462, 0.045),
]
FS8 = ([(n, z, r*f, e*f, "DESI") for n, z, r, e, f in DESI_FS]
       + [(n, z, v, e, "classic") for n, z, v, e in CLASSIC])
FS8.sort(key=lambda t: t[1])

# ----------------------------------------------------------------------------------------------
# 3. Backgrounds: flat LCDM (Om, p) and flat w0waCDM (Om, w0, wa, p), p = c/(H0 r_d)
# ----------------------------------------------------------------------------------------------
def E_of_z(z, Om, w0=-1.0, wa=0.0):
    z = np.asarray(z, float)
    fde = (1+z)**(3*(1+w0+wa)) * np.exp(-3*wa*z/(1+z))     # CPL, smooth DE
    return np.sqrt(Om*(1+z)**3 + (1-Om)*fde)

def bao_chi2(theta, model):
    if model == "lcdm": Om, p = theta; w0, wa = -1.0, 0.0
    else:               Om, w0, wa, p = theta
    if not (0.05 < Om < 0.7): return 1e10
    zg = np.linspace(0, 2.45, 1200)
    invE = 1.0/E_of_z(zg, Om, w0, wa)
    if not np.all(np.isfinite(invE)): return 1e10
    DMg = p*integrate.cumulative_trapezoid(invE, zg, initial=0.0)   # D_M/r_d
    chi2 = 0.0
    for _, z, dv, sig in BAO_DV:
        DM = np.interp(z, zg, DMg); DH = p/E_of_z(z, Om, w0, wa)
        DV = (z*DM**2*DH)**(1/3.)
        chi2 += ((DV-dv)/sig)**2
    for _, z, dm, sm, dh, sh, r in BAO_MH:
        DM = np.interp(z, zg, DMg); DH = p/E_of_z(z, Om, w0, wa)
        d = np.array([DM-dm, DH-dh])
        cov = np.array([[sm**2, r*sm*sh], [r*sm*sh, sh**2]])
        chi2 += d @ np.linalg.solve(cov, d)
    return chi2

def fit_background(model):
    if model == "lcdm":
        bounds = [(0.15, 0.5), (25., 35.)]
    else:  # DESI priors: w0 in [-3,1], wa in [-3,2] (arXiv:2503.14738 Sec. cosmology)
        bounds = [(0.15, 0.5), (-3., 1.), (-3., 2.), (25., 35.)]
    res = optimize.differential_evolution(bao_chi2, bounds, args=(model,), seed=42,
                                          tol=1e-10, maxiter=2000, polish=True)
    return res.x, res.fun

# ----------------------------------------------------------------------------------------------
# 4. Linear growth in GR with smooth (non-clustering) DE:
#    D'' + [2 + dlnE/dlna] D' = 1.5 Om a^-3 E^-2 D    (' = d/dlna)
# ----------------------------------------------------------------------------------------------
def growth(Om, w0=-1.0, wa=0.0):
    lna = np.linspace(np.log(1e-3), 0.0, 4000)
    def rhs(la, y):
        a = np.exp(la); z = 1/a - 1
        E = E_of_z(z, Om, w0, wa)
        dla = 1e-4
        dlnE = (np.log(E_of_z(1/np.exp(la+dla)-1, Om, w0, wa))
                - np.log(E_of_z(1/np.exp(la-dla)-1, Om, w0, wa)))/(2*dla)
        Oma = Om*a**-3/E**2
        return [y[1], -(2+dlnE)*y[1] + 1.5*Oma*y[0]]
    a0 = np.exp(lna[0])
    sol = integrate.solve_ivp(rhs, (lna[0], 0), [a0, a0], t_eval=lna,
                              rtol=1e-8, atol=1e-12, method="RK45")
    D = sol.y[0]; f = sol.y[1]/sol.y[0]
    zs = 1/np.exp(lna) - 1
    return lambda z: (np.interp(np.log(1/(1+np.asarray(z, float))), lna, f)
                      * np.interp(np.log(1/(1+np.asarray(z, float))), lna, D) / D[-1])
    # returns f(z) * D(z)/D(0);  fsigma8_pred = sigma8_0 * that

def null_test(label, theta, model):
    if model == "lcdm": Om, p = theta; w0, wa = -1.0, 0.0
    else:               Om, w0, wa, p = theta
    fD = growth(Om, w0, wa)
    z   = np.array([t[1] for t in FS8]); obs = np.array([t[2] for t in FS8])
    err = np.array([t[3] for t in FS8]); src = np.array([t[4] for t in FS8])
    g = fD(z)                                   # f(z) D(z)/D(0)
    s8 = np.sum(obs*g/err**2)/np.sum(g**2/err**2)   # analytic LSQ amplitude
    pred = s8*g
    O = obs/pred; Oe = err/pred
    print(f"\n--- {label} ---")
    if model == "lcdm":
        print(f"  best fit: Om={Om:.4f}, c/(H0 r_d)={p:.3f}   sigma8_0(fit)={s8:.4f}")
    else:
        print(f"  best fit: Om={Om:.4f}, w0={w0:+.3f}, wa={wa:+.3f}, c/(H0 r_d)={p:.3f}   "
              f"sigma8_0(fit)={s8:.4f}")
    print(f"  {'point':<12}{'z':>7}{'fs8 obs':>10}{'fs8 pred':>10}{'O':>8}{'+/-':>8}{'pull':>7}")
    for i in range(len(z)):
        print(f"  {FS8[i][0]:<12}{z[i]:>7.3f}{obs[i]:>10.3f}{pred[i]:>10.3f}"
              f"{O[i]:>8.3f}{Oe[i]:>8.3f}{(O[i]-1)/Oe[i]:>+7.2f}")
    w = 1/Oe**2
    def band(mask, name):
        if mask.sum() == 0: return
        m = np.sum(O[mask]*w[mask])/np.sum(w[mask]); e = np.sum(w[mask])**-0.5
        print(f"    {name:<38} <O> = {m:.4f} +/- {e:.4f}   ({(m-1)/e:+.2f} sigma from 1, "
              f"N={mask.sum()})")
        return (m-1)/e
    print(f"  weighted mean O (all)          = {np.sum(O*w)/np.sum(w):.4f} "
          f"+/- {np.sum(w)**-0.5:.4f}   [~1 BY CONSTRUCTION: amplitude was fit]")
    chi2 = np.sum(((O-1)/Oe)**2); ndof = len(z)-1
    print(f"  chi2(O=1) = {chi2:.2f} / {ndof} dof   (PTE = {1-stats.chi2.cdf(chi2, ndof):.3f})")
    sig_band = band((z >= 0.3) & (z <= 0.8), "HINT BAND z in [0.3, 0.8]")
    band(z < 0.3,  "z < 0.3")
    band(z > 0.8,  "z > 0.8")
    band((z >= 0.3) & (z <= 0.8) & (src == "DESI"),    "  band, DESI-DR1-FS points only")
    band((z >= 0.3) & (z <= 0.8) & (src == "classic"), "  band, classic points only")
    return chi2, ndof, sig_band

# ================================================================================================
print("="*96)
print("STEP 1 -- fit backgrounds to DESI DR2 BAO alone (13 data points, 7 tracers)")
print("="*96)
th_l, c2_l = fit_background("lcdm")
th_w, c2_w = fit_background("w0wa")
print(f"  flat LCDM  : Om={th_l[0]:.4f}, c/(H0 r_d)={th_l[1]:.3f}"
      f"            chi2={c2_l:.2f} / {13-2} dof")
print(f"  w0waCDM    : Om={th_w[0]:.4f}, w0={th_w[1]:+.3f}, wa={th_w[2]:+.3f}, "
      f"c/(H0 r_d)={th_w[3]:.3f}   chi2={c2_w:.2f} / {13-4} dof")
print(f"  Delta chi2 (LCDM - w0wa) = {c2_l-c2_w:.2f} for 2 extra params")
ok = (th_w[1] > -1) and (th_w[2] < 0) and (0.28 < th_l[0] < 0.32)
print(f"  SANITY GATE (w0>-1, wa<0, LCDM Om in [0.28,0.32]): {'PASS' if ok else 'FAIL -- verdict void'}")

print()
print("="*96)
print("STEP 2 -- growth null test O(z) = fsigma8_obs / fsigma8_pred  (13 points, amplitude free)")
print("="*96)
c2l, ndl, bl = null_test("background: flat LCDM (from BAO alone)", th_l, "lcdm")
c2w, ndw, bw = null_test("background: w0waCDM (from BAO alone)",   th_w, "w0wa")

print()
print("="*96)
print("READ (pre-registered in docstring)")
print("="*96)
verdict_null = (abs(bl) < 2) and (abs(bw) < 2) and (1-stats.chi2.cdf(c2w, ndw) > 0.05)
print(f"  * hint-band z=[0.3,0.8] deviation: {bl:+.2f} sigma (LCDM bg), {bw:+.2f} sigma (w0wa bg)")
print(f"  * chi2(O=1): {c2l:.1f}/{ndl} (LCDM bg), {c2w:.1f}/{ndw} (w0wa bg)")
if verdict_null:
    print("  * VERDICT: NULL -- O(z)=1 within 2 sigma; no perturbation-level support for a COUPLED")
    print("    dark sector. The DESI evolving-DE hint stays background-only, where coupled and")
    print("    uncoupled models are exactly degenerate (Kunz astro-ph/0702615).")
else:
    print("  * VERDICT: DEVIATION -- a coherent growth-vs-background mismatch survives; check")
    print("    subset consistency lines above before believing it. If DESI-only and classic-only")
    print("    agree in sign, this warrants the full paper.")
print("  * KNOWN CAVEATS (cannot be fixed at prototype level):")
print("    - every fsigma8 point was measured against a FIDUCIAL-COSMOLOGY TEMPLATE (AP +")
print("      dilation); the DESI fsig_s8->fsig8 conversion additionally assumes the fiducial r_d.")
print("      A full analysis must redo the compression consistently -- a prototype cannot.")
print("    - cross-bin covariances ignored; DESI DR1 FS and eBOSS overlap in sky/volume, so the")
print("      combined-sample errors are mildly optimistic; asymmetric errors symmetrized, MAP used.")
print("    - amplitude sigma8_0 is marginalized, so this is a SHAPE-only null test; a coupling")
print("      that only rescales sigma8 uniformly in z is invisible here (it is degenerate with")
print("      sigma8_0 anyway without a CMB anchor).")
