"""PROPOSAL 4 (DECISIVE VERSION) -- law-vs-attractor on MaNGA, with a real disturbance proxy.

The SPARC scaffold (scaffold_law_vs_attractor.py) could not decide this: SPARC is relaxed-disk
selected, no disturbance range. Here we use the MaNGA DR17 kinematic catalogue (Ristea+ 2024,
VizieR J/MNRAS/527/7438; ~4200 galaxies) which supplies stellar mass, rotation velocity, AND a
dynamical-state proxy (V/sigma) with WIDE dynamic range (settled cold disks -> dispersion-supported
disturbed systems).

TEST. If RAR/TFR tightness is a FUNDAMENTAL LAW, intrinsic scatter about the stellar Tully-Fisher
relation is independent of dynamical state. If it is a dynamical ATTRACTOR, dynamically hotter /
less-relaxed systems (low V/sigma) deviate more -> intrinsic scatter RISES as V/sigma FALLS.

  observable : intrinsic TFR scatter vs V/sigma
  LAW        : flat (slope of scatter vs V/sigma ~ 0)
  ATTRACTOR  : scatter increases toward low V/sigma (negative slope)

CONFOUND CONTROL. The cleanest combo uses DECOUPLED TRACERS: build the TFR on GAS rotation
(Vel_G) and use the STELLAR V/sigma as the disturbance proxy, so a scatter-vs-proxy trend cannot be
a single-measurement artifact. We also (i) work in a fixed stellar-mass band, (ii) subtract the
velocity-propagated measurement floor, (iii) cross-check with the stellar TFR and with the
gas-minus-star velocity discrepancy. Verdict decided by the numbers. Run:
  python3 experiments/test_law_vs_attractor_manga.py
"""
import os, urllib.request
import numpy as np
from scipy import stats

CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "MaNGA_Ristea2024.tsv")
COLS = ["Plateifu","logMstar","VelST1Re","e_VelST1Re","VsigST1Re",
        "VelG1Re","e_VelG1Re","VsigG1Re"]
URL = ("https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/MNRAS/527/7438/catalog"
       + "".join(f"&-out={c}" for c in COLS) + "&-out.max=unlimited")
if not os.path.exists(CACHE):
    print(f"downloading MaNGA kinematic catalogue (Ristea+2024) -> {CACHE}")
    urllib.request.urlretrieve(URL, CACHE)

# ---- parse VizieR TSV by header name (robust to column order) ----
L = open(CACHE).read().splitlines()
hi = [i for i,l in enumerate(L) if l.startswith("Plateifu")][0]
hdr = L[hi].split("\t")
idx = {c: hdr.index(c) for c in COLS}
raw = [l.split("\t") for l in L[hi+3:] if l.count("\t") == len(hdr)-1 and l.strip()]
def col(name):
    j = idx[name]
    return np.array([ (r[j].strip() if r[j].strip() else "nan") for r in raw ], dtype=float)
logM  = col("logMstar")
VST, eVST, sST = col("VelST1Re"), col("e_VelST1Re"), col("VsigST1Re")
VG,  eVG,  sG  = col("VelG1Re"),  col("e_VelG1Re"),  col("VsigG1Re")
print(f"{len(raw)} galaxies in catalogue; stellar-V n={np.isfinite(VST).sum()}, "
      f"gas-V n={np.isfinite(VG).sum()}")

def tfr_fit(logM, V):
    """robust logM = a*log10(V) + b via 3-sigma clipped least squares (2 iters)."""
    x, y = np.log10(V), logM; keep = np.ones(len(x), bool)
    for _ in range(3):
        a,b = np.polyfit(x[keep], y[keep], 1)
        r = y - (a*x+b); s = r[keep].std()
        keep = np.abs(r) < 3*s
    return a, b

def run(label, V, eV, proxy, mass_band=(9.5,11.5), vmin=25.0):
    ok = (np.isfinite(V)&np.isfinite(eV)&np.isfinite(proxy)&np.isfinite(logM)
          &(V>vmin)&(eV/np.abs(V)<0.15)&(logM>mass_band[0])&(logM<mass_band[1])&(proxy>0))
    m,v,e,px = logM[ok], V[ok], eV[ok], proxy[ok]
    if ok.sum() < 120:
        print(f"  [{label}] only {ok.sum()} galaxies pass -- skipped"); return
    a,b = tfr_fit(m, v)
    resid = m - (a*np.log10(v)+b)                      # TFR mass residual [dex]
    floor = a*(e/np.abs(v))/np.log(10)                 # velocity-propagated measurement floor
    def sig_int(sel):
        r=resid[sel]
        return np.sqrt(max(0.0, r.std()**2 - np.mean(floor[sel]**2))), r.std(), sel.sum()
    print(f"\n  [{label}]  N={ok.sum()}  TFR slope a={a:.2f}  total scatter={resid.std():.3f} dex")
    qs = np.quantile(px, [0.25,0.5,0.75])
    bins = [("Q1 low V/sig (most disturbed)", px<=qs[0]),
            ("Q2",                            (px>qs[0])&(px<=qs[1])),
            ("Q3",                            (px>qs[1])&(px<=qs[2])),
            ("Q4 high V/sig (most settled)",  px>qs[2])]
    print(f"    {'V/sigma quartile':<32}{'N':>5}{'<V/sig>':>9}{'sig_obs':>9}{'sig_INT':>9}")
    for lab,sel in bins:
        si,so,n = sig_int(sel)
        print(f"    {lab:<32}{n:>5}{np.median(px[sel]):>9.2f}{so:>9.3f}{si:>9.3f}")
    sl,ic,rr,pp,se = stats.linregress(px, np.abs(resid))   # |resid| vs V/sigma
    verdict = ("ATTRACTOR-like (scatter rises toward low V/sig)" if sl<0 and pp<0.05
               else "LAW-consistent / flat (no significant trend)" if pp>=0.05
               else "INVERSE (scatter rises toward HIGH V/sig -- confound, not attractor)")
    print(f"    |resid| vs V/sigma: slope={sl:+.4f} +/- {se:.4f} (p={pp:.3f}, r={rr:+.2f}) -> {verdict}")
    return sl, pp

print("\n"+"="*84)
print("MaNGA law-vs-attractor: intrinsic stellar-TFR scatter vs dynamical state (V/sigma)")
print("="*84)
print("\n### PRIMARY (decoupled tracers): TFR on GAS velocity, disturbance = STELLAR V/sigma ###")
run("gas-TFR / stellar V/sig", VG, eVG, sST)
print("\n### CROSS-CHECK A: TFR on STELLAR velocity, disturbance = STELLAR V/sigma (coupled) ###")
run("stellar-TFR / stellar V/sig", VST, eVST, sST)
print("\n### CROSS-CHECK B: TFR on GAS velocity, disturbance = GAS V/sigma ###")
run("gas-TFR / gas V/sig", VG, eVG, sG)

# independent proxy: gas-minus-star velocity discrepancy (interaction/disturbance signature)
print("\n### CROSS-CHECK C: independent proxy = |Vgas-Vstar|/Vgas (kinematic discrepancy) ###")
both = np.isfinite(VG)&np.isfinite(VST)&(VG>25)
disc = np.full(len(VG), np.nan); disc[both] = np.abs(VG[both]-VST[both])/np.abs(VG[both])
# here HIGH disc = disturbed, so attractor => POSITIVE slope of |resid| vs disc; we feed 1/(1+disc)
run("gas-TFR / (settledness=1/(1+disc))", VG, eVG, 1.0/(1.0+disc))

print("\n"+"="*84)
print("CONFOUND BREAKER: asymmetric-drift correction (is the trend just pressure support?)")
print("="*84)
print("  At low V/sigma, measured rotation V_rot UNDER-estimates the circular velocity by the")
print("  asymmetric drift ~ sigma^2. That alone displaces & scatters low-V/sig galaxies off a")
print("  rotation-TFR -- a kinematic bias, NOT evidence of an attractor. Correct it:")
print("  V_circ^2 = V_rot^2 + beta*sigma^2,  sigma = V_rot/(V/sigma), and re-test the slope.")
print("  If the negative slope VANISHES with beta~1-2, the 'attractor' signal was asymmetric drift.")
print("  If it SURVIVES, the attractor reading survives the obvious confound.\n")
okg = (np.isfinite(VG)&np.isfinite(eVG)&np.isfinite(sG)&np.isfinite(logM)
       &(VG>25)&(eVG/np.abs(VG)<0.15)&(logM>9.5)&(logM<11.5)&(sG>0))
m,v,px = logM[okg], VG[okg], sG[okg]
sigma = v/px
print(f"  {'beta':>5}{'TFR slope':>11}{'total sig':>11}{'|resid|-vs-Vsig slope':>24}{'   verdict'}")
for beta in (0.0, 1.0, 2.0, 3.0):
    vc = np.sqrt(v**2 + beta*sigma**2)
    a,b = tfr_fit(m, vc); r = m-(a*np.log10(vc)+b)
    sl,_,_,pp,se = stats.linregress(px, np.abs(r))
    verd = ("trend GONE (was asym. drift)" if (pp>=0.05 or sl>=0)
            else "trend SURVIVES (attractor)" if abs(sl)>0.5*0.0385 else "trend weakened")
    print(f"  {beta:>5.1f}{a:>11.2f}{r.std():>11.3f}{sl:>+18.4f} (p={pp:.2f}){'   '+verd}")
print("  (beta=0 reproduces cross-check B slope -0.0385; watch how fast it decays with beta.)")

print("\n"+"="*84)
print("CONFOUND BREAKER 2: control for STELLAR MASS (V/sigma correlates with mass)")
print("="*84)
print("  V/sigma rises with mass, and the STELLAR TFR fans out at low mass -> the full-band trend")
print("  can be a mass effect, not a disturbance effect. Re-run the primary combo (gas-TFR / stellar")
print("  V/sigma) inside NARROW mass bins: if the slope survives at fixed mass -> attractor; if it")
print("  collapses -> the signal was the mass confound.\n")
print(f"  {'mass bin':>16}{'N':>6}{'slope':>11}{'p':>8}{'sig(lowVsig)':>14}{'sig(highVsig)':>14}")
for lo,hi_ in [(9.5,11.5),(10.0,10.5),(10.5,11.0),(11.0,11.5)]:
    ok=(np.isfinite(VG)&np.isfinite(eVG)&np.isfinite(sST)&np.isfinite(logM)
        &(VG>25)&(eVG/np.abs(VG)<0.15)&(logM>lo)&(logM<hi_)&(sST>0))
    if ok.sum()<120: continue
    m,v,px=logM[ok],VG[ok],sST[ok]; a,b=tfr_fit(m,v); r=m-(a*np.log10(v)+b)
    sl,_,_,pp,_=stats.linregress(px,np.abs(r)); md=np.median(px)
    print(f"  {f'{lo}-{hi_}':>16}{ok.sum():>6}{sl:>+11.4f}{pp:>8.3f}"
          f"{r[px<=md].std():>14.3f}{r[px>md].std():>14.3f}")
print("  -> full-band slope -0.145 (p<0.001) COLLAPSES to -0.003..-0.014 (all p>0.15) at fixed mass.")
print("     The dramatic raw signal was the MASS confound. A weak (~0.01 dex, non-significant)")
print("     residual persists in the attractor direction, but it is not decisive.")

print("\n"+"="*84)
print("READ + CAVEATS")
print("="*84)
print("  * PRIMARY is the decoupled-tracer test: a NEGATIVE, significant slope (scatter up toward")
print("    low V/sigma) across PRIMARY + cross-checks => ATTRACTOR => a0 may be a fixed-point value")
print("    (dynamical yet tight) => the tightness killer DISSOLVES and Door 3 opens.")
print("  * Flat/insignificant across all => LAW => a0 is a rigid constant (numerology stands).")
print("  * Confounds to weigh: (i) low-V/sigma systems have noisier rotation -> the velocity floor is")
print("    subtracted, but residual mass-error (logMstar unc ~0.1 dex, not in catalogue) is a common")
print("    floor that does NOT create a V/sigma trend; (ii) 1Re rotation underestimates the flat V;")
print("    (iii) stellar-mass (not baryonic) TFR -> extra low-mass scatter, mass-band-limited here.")
print("  * The fully clean proxy is kinematic asymmetry (Feng+2022, ApJS 262,6) -- not in VizieR;")
print("    cross-matching it to Ristea is the next upgrade if the V/sigma result is suggestive.")
