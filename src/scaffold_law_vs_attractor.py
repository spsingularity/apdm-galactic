"""PROPOSAL 4 (SCAFFOLD) -- is RAR tightness a fundamental LAW or a dynamical ATTRACTOR?

This is the move that decides whether the a0 dilemma is solvable. Two interpretations of the
RAR's tightness make OPPOSITE predictions for how intrinsic scatter depends on dynamical state:

  * FUNDAMENTAL LAW (MOND-as-axiom): the relation holds always, even out of equilibrium.
    Intrinsic scatter is pure noise -> INDEPENDENT of relaxation state. => a0 must be a rigid
    constant => stuck with numerology (Proposal 1's fate).
  * ATTRACTOR (RG fixed point / dissipative self-organization): the relation is dynamically
    enforced; disturbed systems deviate and relax back. Intrinsic scatter DECREASES with
    relaxation (more completed orbits / longer since disturbance). => a0 can be the fixed-point
    value -- DYNAMICAL WITHOUT SCATTER. The tightness killer dissolves.

  DISCRIMINATING PREDICTION
  -------------------------------------------------------------------------------------
  observable: intrinsic RAR scatter vs dynamical-relaxation proxy, at fixed acceleration regime
  LAW      -> flat (slope 0 within errors)
  ATTRACTOR-> decreasing (scatter smaller for more-relaxed systems)
  -------------------------------------------------------------------------------------

WHAT THIS SCAFFOLD DOES (in-hand data): a FIRST LOOK using the only relaxation proxy available in
SPARC -- the number of completed orbits N_orb = t_age / t_dyn(r), t_dyn=2*pi*r/V(r). Inner/faster
points have relaxed longer. We control for the acceleration regime (the obvious confound) by
working inside the deep-MOND band and inside a narrow g_bar slice.

WHAT IT CANNOT DO: SPARC is a RELAXED-DISK-selected sample -> tiny dynamic range in disturbance,
and N_orb is confounded with radius/regime. The DECISIVE version needs disturbed systems (see the
DATA GAP block). This scaffold sets up the machinery, gives the SPARC first look, and specifies
the external data. Run: python3 experiments/scaffold_law_vs_attractor.py
"""
import os, urllib.request
import numpy as np
from scipy import optimize, stats

CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "SPARC_MassModels.mrt")
URL = "https://raw.githubusercontent.com/jobovy/sparc-rotation-curves/main/data/MassModels_Lelli2016c.mrt"
if not os.path.exists(CACHE):
    urllib.request.urlretrieve(URL, CACHE)

YD, YB = 0.5, 0.7
CONV = (1e3)**2/3.0856775814913673e19            # (km/s)^2/kpc -> m/s^2
KPC_KMS_TO_GYR = 3.0856775814913673e16/3.1557e16 # (kpc/(km/s)) -> Gyr
T_AGE = 10.0                                       # Gyr (proxy for available relaxation time)

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
R,Vobs,eV,gbar,gobs = R[mask],Vobs[mask],eV[mask],gbar[mask],gobs[mask]

rar=lambda gb,gd: gb/(1.0-np.exp(-np.sqrt(gb/gd)))
res=optimize.minimize_scalar(lambda lg: np.sum((np.log10(rar(gbar,10**lg))-np.log10(gobs))**2),
                             bounds=(-10.5,-9.0),method='bounded')
gd=10**res.x
resid = np.log10(gobs) - np.log10(rar(gbar,gd))    # RAR residual [dex]
sig_meas = 2.0*(eV/Vobs)/np.log(10)                # per-point measurement floor in log gobs [dex]
N_orb = T_AGE/(2*np.pi*R/Vobs*KPC_KMS_TO_GYR)      # completed orbits (relaxation proxy)

def intrinsic(sel):
    r=resid[sel]; m=sig_meas[sel]
    if len(r)<15: return np.nan,np.nan,len(r)
    return np.sqrt(max(0.0, r.std()**2 - np.mean(m**2))), r.std(), len(r)

print("="*84)
print("PROPOSAL 4 SCAFFOLD: law-vs-attractor via intrinsic RAR scatter vs N_orb (relaxation proxy)")
print("="*84)
print(f"  g_dagger={gd:.3e} m/s^2; {mask.sum()} pts; median N_orb={np.median(N_orb):.0f}, "
      f"range [{np.percentile(N_orb,5):.0f}, {np.percentile(N_orb,95):.0f}]")

print("\n-- (1) deep-MOND band (g_bar < 0.1 g_dagger), split by N_orb tertiles --")
dm = gbar < 0.1*gd
q = np.quantile(N_orb[dm], [1/3, 2/3])
print(f"  {'N_orb bin':<22}{'N':>6}{'sigma_obs':>11}{'sigma_meas':>12}{'sigma_INT':>11}")
for lab, sel in [("least relaxed (low)", dm&(N_orb<=q[0])),
                 ("middle",             dm&(N_orb>q[0])&(N_orb<=q[1])),
                 ("most relaxed (high)",dm&(N_orb>q[1]))]:
    si,so,n = intrinsic(sel)
    mm = np.mean(sig_meas[sel])
    print(f"  {lab:<22}{n:>6}{so:>11.3f}{mm:>12.3f}{si:>11.3f}")

print("\n-- (2) confound control: narrow g_bar slice (0.02-0.06 g_dagger), scatter vs N_orb --")
sl = dm & (gbar>0.02*gd) & (gbar<0.06*gd)
if sl.sum()>=45:
    qq=np.quantile(N_orb[sl],[1/3,2/3])
    for lab,sel in [("low N_orb",sl&(N_orb<=qq[0])),("mid",sl&(N_orb>qq[0])&(N_orb<=qq[1])),
                    ("high N_orb",sl&(N_orb>qq[1]))]:
        si,so,n=intrinsic(sel)
        print(f"  {lab:<12}{n:>6} pts  sigma_obs={so:.3f}  sigma_INT={si:.3f}")
    sl_slope,_,sl_r,sl_p,_ = stats.linregress(np.log10(N_orb[sl]), np.abs(resid[sl]))
    print(f"  |residual| vs log N_orb (in slice): slope={sl_slope:+.3f} (p={sl_p:.2f}) "
          f"-> {'DECREASING (attractor-like)' if sl_slope<0 and sl_p<0.05 else 'flat/inconclusive (law-consistent or underpowered)'}")

print("\n" + "="*84)
print("FIRST-LOOK READ + WHY IT IS NOT DECISIVE")
print("="*84)
print("  * Interpret slope of sigma_INT (or |resid|) vs N_orb: <0 & significant => attractor;")
print("    ~0 => law-consistent OR underpowered. SPARC almost certainly gives the latter because:")
print("    (a) it is relaxed-disk-SELECTED -> little disturbance range; (b) N_orb correlates with")
print("    radius/regime; (c) intrinsic scatter here is already near the measurement floor.")
print("  * So treat this as machinery validation + a null, NOT a verdict.")
print()
print("="*84)
print("DATA GAP -- what the DECISIVE test needs (all existing or near-term, no DM detection)")
print("="*84)
print("  relaxation proxy with real dynamic range, then measure intrinsic RAR scatter vs proxy:")
print("   1. MaNGA/SAMI IFU kinematic asymmetry (v_asym, kinemetry) -> disturbed vs settled disks.")
print("   2. Post-merger catalogs (time-since-merger) -> scatter vs merger age; attractor predicts")
print("      elevated scatter decaying over a relaxation time (~Gyr), law predicts none.")
print("   3. Tidal dwarf galaxies (dynamically young, DM-poor) -> a clean out-of-equilibrium probe:")
print("      LAW says they sit exactly on the RAR; ATTRACTOR says they scatter high and relax on.")
print("  Decision rule: intrinsic scatter falls with relaxation across >=2 of these => ATTRACTOR =>")
print("  a0 may be the fixed-point value (dynamical, tight) => the tightness killer is DISSOLVED, and")
print("  Door 3 (protected-structure a0) opens. Flat across all => LAW => a0 is a constant (numerology")
print("  stands) and no mechanism will explain its value -- stop building a0-generating mechanisms.")
