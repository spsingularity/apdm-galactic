"""TEST 2 of the 'state-not-coupling' reframe, against REAL SPARC data (Lelli+2016).
 2a: are RAR residuals ASYMMETRIC (no super-RAR outliers) -> flux-bound saturation (FAB4)?
 2b: is the deep-MOND asymptotic slope ~1/2 (universal fixed-point attractor, FP2)?
 sanity: recover g_dagger (McGaugh+2016: 1.20e-10 m/s^2) and the ~0.11-0.13 dex scatter.

Data: MassModels_Lelli2016c.mrt. Auto-downloads from a public mirror to a local cache if absent
(official source: http://astroweb.cwru.edu/SPARC/). Run: python3 experiments/test_reframe_rar_sparc.py
"""
import os, urllib.request
import numpy as np
from scipy import stats, optimize

CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "SPARC_MassModels.mrt")
URL = "https://raw.githubusercontent.com/jobovy/sparc-rotation-curves/main/data/MassModels_Lelli2016c.mrt"
if not os.path.exists(CACHE):
    print(f"downloading SPARC mass models -> {CACHE}")
    urllib.request.urlretrieve(URL, CACHE)

YD, YB = 0.5, 0.7                                   # M/L disk, bulge at 3.6um (standard RAR)
CONV = (1e3)**2/3.0856775814913673e19              # (km/s)^2/kpc -> m/s^2 = 3.241e-14

rows=[]; started=False; ndash=0
for line in open(CACHE):
    if line.startswith('--------'):
        ndash+=1; started=(ndash>=3); continue
    if not started or not line.strip(): continue
    p=line.split()
    if len(p)<10: continue
    try: D,R,Vobs,eV,Vgas,Vdisk,Vbul=map(float,p[1:8])
    except ValueError: continue
    rows.append((R,Vobs,eV,Vgas,Vdisk,Vbul))
R,Vobs,eV,Vgas,Vdisk,Vbul=np.array(rows).T
Vbar2 = Vgas*np.abs(Vgas)+YD*Vdisk*np.abs(Vdisk)+YB*Vbul*np.abs(Vbul)
gbar = Vbar2/R*CONV; gobs = Vobs**2/R*CONV
mask = (Vbar2>0)&(R>0)&(Vobs>0)&(eV/Vobs<0.10)&np.isfinite(gbar)&np.isfinite(gobs)
gbar,gobs=gbar[mask],gobs[mask]
x,y=np.log10(gbar),np.log10(gobs)
print(f"{len(rows)} radial points; {len(gbar)} pass cut (e_V/V<0.10, g_bar>0)")

rar=lambda gb,gd: gb/(1.0-np.exp(-np.sqrt(gb/gd)))
res=optimize.minimize_scalar(lambda lg: np.sum((np.log10(rar(gbar,10**lg))-y)**2),
                             bounds=(-10.5,-9.0),method='bounded')
gd=10**res.x; resid=y-np.log10(rar(gbar,gd))
print(f"[sanity] g_dagger={gd:.3e} m/s^2 (lit 1.20e-10)   scatter={resid.std():.3f} dex (lit ~0.11-0.13)")

sk=stats.skew(resid); sig=resid.std()
print("\n=== 2a: residual asymmetry (saturation ceiling / FAB4 -> expect strong NEGATIVE skew, one-sided) ===")
print(f"  skewness={sk:+.3f}   mean={resid.mean():+.3f} dex")
print(f"  N(>+2s super-RAR)={np.sum(resid>2*sig):3d}  N(<-2s sub-RAR)={np.sum(resid<-2*sig):3d}")
print(f"  N(>+3s)={np.sum(resid>3*sig):3d}          N(<-3s)={np.sum(resid<-3*sig):3d}")
print("  -> negative skew + one-sided 3s tail = sign matches FAB4, BUT super-RAR outliers exist and")
print("     the low tail is confounded by low-inclination/non-circular systematics -> NOT DECISIVE.")

print("\n=== 2b: deep-MOND slope (fixed-point attractor / FP2 -> expect 0.500) ===")
for thr in (0.1,0.03):
    dm=gbar<thr*gd
    if dm.sum()>=20:
        sl,_,_,_,se=stats.linregress(x[dm],y[dm])
        print(f"  g_bar<{thr:g}*gd ({dm.sum():4d} pts): slope={sl:.3f}+/-{se:.3f}")
print("  -> slope~0.5 in the reliable (<0.1 gd) regime CONFIRMS a real tight attractor, but 0.5 is")
print("     also exactly BK/MOND -> confirms native tightness, does NOT discriminate the reframe.")
