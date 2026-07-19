"""PROPOSAL 4, CLEAN PROXY -- law-vs-attractor with an INDEPENDENT disturbance measure.

The MaNGA V/sigma test (test_law_vs_attractor_manga.py) showed a strong attractor-direction signal
that COLLAPSED under stellar-mass control -> the raw trend was a mass confound. The clean follow-up
needs a disturbance proxy that is NOT the galaxy's own internal kinematics. The obvious choice was
Feng+2022 kinematic asymmetry, but that per-galaxy catalog was never publicly released (only summary
tables). SUBSTITUTE: the GEMA-VAC (Argudo-Fernandez+; SDSS DR17), which gives an EXTERNAL disturbance
measure -- the tidal strength of the nearest neighbour Q_nn (log ratio of external tidal force to
internal binding) -- plus a merger probability p_merger for close pairs, all keyed by MaNGA ID and
directly joinable to the Ristea+2024 kinematic catalogue. This is the same external-tidal proxy Feng
himself used (his Q_local).

TEST (unchanged): does intrinsic stellar-TFR scatter rise with disturbance?
  LAW       -> flat at fixed mass;  ATTRACTOR -> scatter rises with tidal strength / merger prob.
CENTERPIECE: control for stellar mass (the confound that killed the V/sigma signal). Q_nn correlates
with mass (massive galaxies sit in denser environments), so the test MUST be done at fixed mass.

Data (auto-download to caches): Ristea+2024 (VizieR J/MNRAS/527/7438); GEMA_2.0.2.fits (SDSS DR17).
Run: python3 experiments/test_law_vs_attractor_gema.py   (needs astropy)
"""
import os, urllib.request
import numpy as np
from scipy import stats
from astropy.io import fits

HERE = os.path.dirname(os.path.abspath(__file__))
RCACHE = os.path.join(HERE, "data", "MaNGA_Ristea2024_id.tsv")
GCACHE = os.path.join(HERE, "data", "GEMA_2.0.2.fits")
RCOLS = ["MaNGA","logMstar","VelG1Re","e_VelG1Re","VsigST1Re","VelST1Re","e_VelST1Re"]
RURL = ("https://vizier.cds.unistra.fr/viz-bin/asu-tsv?-source=J/MNRAS/527/7438/catalog"
        + "".join(f"&-out={c}" for c in RCOLS) + "&-out.max=unlimited")
GURL = "https://data.sdss.org/sas/dr17/manga/gema/2.0.2/GEMA_2.0.2.fits"
if not os.path.exists(RCACHE): print("downloading Ristea+2024 ..."); urllib.request.urlretrieve(RURL, RCACHE)
if not os.path.exists(GCACHE): print("downloading GEMA-VAC ...");     urllib.request.urlretrieve(GURL, GCACHE)

# ---- Ristea ----
L = open(RCACHE).read().splitlines()
hi = [i for i,l in enumerate(L) if l.startswith("MaNGA")][0]; hdr = L[hi].split("\t")
ix = {c: hdr.index(c) for c in hdr}
raw = [l.split("\t") for l in L[hi+3:] if l.count("\t") == len(hdr)-1 and l.strip()]
def rc(n):
    j = ix[n]; return np.array([(r[j].strip() or "nan") for r in raw])
mid  = np.array([r[ix["MaNGA"]].strip() for r in raw])
logM = rc("logMstar").astype(float)
VG, eVG = rc("VelG1Re").astype(float), rc("e_VelG1Re").astype(float)

# ---- GEMA: tidal strength (HDU2, 1Mpc all) and merger pairs (HDU12) ----
g = fits.open(GCACHE)
gm = np.char.strip(g[2].data["mangaid"].astype(str))
Qnn_map = dict(zip(gm, np.array(g[2].data["Q_nn"], float)))
Q = np.array([Qnn_map.get(m, np.nan) for m in mid])          # tidal strength of nearest neighbour
pm = np.char.strip(g[12].data["mangaid"].astype(str))
pmerg_map = dict(zip(pm, np.array(g[12].data["p_merger"], float)))
pmerg = np.array([pmerg_map.get(m, np.nan) for m in mid])    # merger prob (close pairs only)
print(f"Ristea {len(mid)} gal; matched tidal Q_nn: {np.isfinite(Q).sum()}; in close-pair table: "
      f"{np.isfinite(pmerg).sum()}")

def tfr_fit(m, v):
    x, y = np.log10(v), m; k = np.ones(len(x), bool)
    for _ in range(3):
        a,b = np.polyfit(x[k], y[k], 1); r = y-(a*x+b); k = np.abs(r) < 3*r[k].std()
    return a, b

base = (np.isfinite(VG)&np.isfinite(eVG)&np.isfinite(logM)&(VG>25)&(eVG/np.abs(VG)<0.15)
        &(logM>9.5)&(logM<11.5))
a,b = tfr_fit(logM[base], VG[base])
resid = logM - (a*np.log10(VG)+b)
floor = a*(eVG/np.abs(VG))/np.log(10)
def sig_int(sel):
    r=resid[sel]; return np.sqrt(max(0.0, r.std()**2 - np.mean(floor[sel]**2))), sel.sum()

print("\n"+"="*82)
print("TEST 1 -- tidal strength Q_nn (external disturbance).  Attractor: scatter UP with Q_nn.")
print("="*82)
ok = base & np.isfinite(Q) & (Q>-90)
print(f"  gas-TFR slope a={a:.2f}; N with valid Q_nn = {ok.sum()}")
print("\n  (1a) NO mass control -- full sample, Q_nn tertiles:")
qs = np.quantile(Q[ok], [1/3,2/3])
for lab,sel in [("low Q_nn (isolated)", ok&(Q<=qs[0])),("mid",ok&(Q>qs[0])&(Q<=qs[1])),
                ("high Q_nn (disturbed)",ok&(Q>qs[1]))]:
    si,n=sig_int(sel); print(f"    {lab:<24} N={n:4d}  <Q_nn>={np.median(Q[sel]):+.2f}  sig_INT={si:.3f}")
sl,_,_,pp,se = stats.linregress(Q[ok], np.abs(resid[ok]))
print(f"    |resid| vs Q_nn slope={sl:+.4f}+/-{se:.4f} (p={pp:.3f})  [uncontrolled]")

print("\n  (1b) MASS-CONTROLLED -- slope of |resid| vs Q_nn inside 0.5-dex mass bins:")
print(f"    {'mass bin':>12}{'N':>6}{'slope':>10}{'p':>8}{'sig(loQ)':>10}{'sig(hiQ)':>10}")
wsum=0; wn=0
for lo in (9.5,10.0,10.5,11.0):
    sel = ok & (logM>=lo) & (logM<lo+0.5)
    if sel.sum()<80: continue
    q=Q[sel]; r=np.abs(resid[sel]); md=np.median(q)
    s,_,_,p,e = stats.linregress(q, r)
    slo,_=sig_int(sel&(Q<=md)); shi,_=sig_int(sel&(Q>md))
    print(f"    {f'{lo}-{lo+0.5}':>12}{sel.sum():>6}{s:>+10.4f}{p:>8.3f}{slo:>10.3f}{shi:>10.3f}")
    wsum += s*sel.sum(); wn += sel.sum()
print(f"    mass-stacked mean slope (N-weighted) = {wsum/wn:+.4f}")
print(f"    (attractor => POSITIVE slope: scatter rises toward HIGH Q_nn; ~0 & p>0.05 => LAW)")

print("\n"+"="*82)
print("TEST 2 -- merger probability p_merger (GEMA close pairs).  Attractor: scatter UP with p_merger.")
print("="*82)
# at fixed mass (10-11.3, best populated), compare close-pair galaxies vs isolated
mb = base & (logM>10.0) & (logM<11.3)
inpair = mb & np.isfinite(pmerg)                                  # in GEMA close-pair table
iso    = mb & ~np.isfinite(pmerg) & np.isfinite(Q) & (Q<np.nanmedian(Q[ok]))  # low-tidal, no pair
for lab,sel in [("close-pair galaxies (disturbed)", inpair),
                ("  of which p_merger>=0.5",          inpair & (pmerg>=0.5)),
                ("isolated (low tidal, no pair)",     iso)]:
    si,n=sig_int(sel)
    print(f"    {lab:<34} N={n:4d}  sig_INT={si:.3f}")
print("    (indicative: close-pair sig_INT vs isolated at same mass; high-p_merger bin tiny after cuts)")

print("\n"+"="*82)
print("READ")
print("="*82)
print("  * Q_nn is an EXTERNAL, kinematics-independent disturbance proxy (unlike V/sigma), and it is")
print("    the tidal measure Feng used. It correlates with mass (r~0.24), so 1b (mass-controlled) is")
print("    the verdict line, not 1a.")
print("  * mass-stacked slope ~0 & bin p-values >0.05 => LAW (tightness fundamental; a0 a constant).")
print("    consistently POSITIVE, significant => ATTRACTOR (a0 may be a dynamical fixed point).")
print("  * Test 2 cross-checks with an orthogonal (merger-pair) proxy; treat as indicative (small N).")
print("  * Caveats carried over: 1Re rotation, stellar-mass (not baryonic) TFR, logMstar unc ~0.1 dex")
print("    as a common floor (does not create a within-mass Q_nn trend).")
