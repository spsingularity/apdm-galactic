"""PROPOSAL 2A prototype -- a METALLICITY INSTRUMENTAL-VARIABLE test of the Gaia wide-binary
gravitational anomaly (Chae 2023, arXiv:2309.10404, claims a ~15-20% velocity boost at
accelerations < 1e-9 m/s^2; Banik et al. find Newtonian; the unresolved crux is hidden-companion
contamination).

THE IV LOGIC. Stellar multiplicity -- and therefore hidden-tertiary contamination -- rises strongly
toward LOW metallicity (the close-companion fraction roughly doubles from [Fe/H]=0 to -1;
Moe, Kratter & Badenes 2019, arXiv:1808.02116). GRAVITY cannot depend on [Fe/H]. So bin the
low-acceleration velocity-excess statistic by metallicity:
  - excess RISES toward low [Fe/H] with the multiplicity-predicted slope  => hidden companions;
  - excess FLAT in [Fe/H]                                                 => gravitational (not contamination).
[Fe/H] is an instrument: it moves the confounder (contamination) but not the hypothesis (gravity).

DATA (fully reproducible, no 1.4 GB catalogue download): a wide-binary sample is built directly from
the Gaia DR3 archive via a TAP self-join (nearby, high-quality pairs), with GSP-Phot [M/H] pulled in
the SAME query. Default volume d < 50 pc (parallax 20-120 mas) -- the largest that the PUBLIC Gaia
TAP self-join runs reliably in one sync call (~16 s); deeper volumes exceed the server's self-join
limit (documented in the note). Cache: experiments/wb_gaia_pairs.csv (gitignored).

PRE-REGISTERED READ (written before computing):
  * significant NEGATIVE slope d(v~)/d[Fe/H] in the low-g regime, matching the contamination sign
    => the anomaly is hidden companions;
  * slope consistent with 0 => supports a gravitational origin (or at least: not contamination);
  * |slope| error too large to distinguish => UNDERPOWERED, report the N needed.
SANITY GATES: high-acceleration control must be (i) flat in [Fe/H] and (ii) have median v~ well
below 1 with ~no v~>1 tail (bound Newtonian orbits); if not, the pipeline is biased -- debug first.
Run: python3 experiments/test_wb_metallicity_iv.py
"""
import os, urllib.request, urllib.parse
import numpy as np
from scipy import stats

CACHE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data", "wb_gaia_pairs.csv")
ADQL = """SELECT a.source_id AS id1,b.source_id AS id2,a.ra AS ra1,a.dec AS de1,
 a.parallax AS plx1,b.parallax AS plx2,a.pmra AS pmra1,a.pmdec AS pmde1,b.pmra AS pmra2,b.pmdec AS pmde2,
 a.pmra_error AS epmra1,a.pmdec_error AS epmde1,b.pmra_error AS epmra2,b.pmdec_error AS epmde2,
 a.phot_g_mean_mag AS g1,b.phot_g_mean_mag AS g2,a.mh_gspphot AS mh1,b.mh_gspphot AS mh2,
 DISTANCE(POINT('ICRS',a.ra,a.dec),POINT('ICRS',b.ra,b.dec))*3600 AS sep_asec
FROM gaiadr3.gaia_source AS a JOIN gaiadr3.gaia_source AS b
 ON 1=CONTAINS(POINT('ICRS',a.ra,a.dec),CIRCLE('ICRS',b.ra,b.dec,0.18))
WHERE a.source_id<b.source_id AND a.parallax BETWEEN 20 AND 120 AND b.parallax BETWEEN 20 AND 120
 AND ABS(a.parallax-b.parallax)<1.5 AND a.parallax_over_error>20 AND b.parallax_over_error>20
 AND a.ruwe<1.4 AND b.ruwe<1.4 AND ABS(a.pmra-b.pmra)<6 AND ABS(a.pmdec-b.pmdec)<6"""
if not os.path.exists(CACHE):
    print("querying Gaia DR3 TAP for wide-binary pairs + [M/H] (d<50pc, ~16s) ...")
    body = urllib.parse.urlencode({"REQUEST":"doQuery","LANG":"ADQL","FORMAT":"csv","QUERY":ADQL}).encode()
    req = urllib.request.Request("https://gea.esac.esa.int/tap-server/tap/sync", data=body)
    with open(CACHE,"wb") as f: f.write(urllib.request.urlopen(req).read())

# ---- Pecaut & Mamajek (2013) main-sequence M_G -> mass (Msun) ----
_MG=np.array([2.2,2.8,3.5,4.1,4.7,5.5,6.7,7.5,8.8,9.8,10.9,12.0,13.0,14.0,15.0])
_MS=np.array([1.55,1.4,1.20,1.08,1.00,0.90,0.80,0.70,0.60,0.50,0.42,0.32,0.22,0.16,0.12])
mass_from_MG=lambda MG: np.interp(MG,_MG,_MS,left=1.7,right=0.09)
G=6.674e-11; Msun=1.989e30; AU=1.496e11

d=np.genfromtxt(CACHE,delimiter=',',names=True)
plx=0.5*(d['plx1']+d['plx2']); dist=1000.0/plx
MG1=d['g1']+5+5*np.log10(plx/1000.0); MG2=d['g2']+5+5*np.log10(plx/1000.0)
Mtot=mass_from_MG(MG1)+mass_from_MG(MG2)
s_AU=d['sep_asec']*dist; s=s_AU*AU
dpm=np.hypot(d['pmra1']-d['pmra2'],d['pmde1']-d['pmde2'])
edpm=np.sqrt(d['epmra1']**2+d['epmde1']**2+d['epmra2']**2+d['epmde2']**2)
dv=4.74047*dpm*(dist/1000.0)*1e3; edv=4.74047*edpm*(dist/1000.0)*1e3    # m/s
vc=np.sqrt(G*Mtot*Msun/s)                                              # m/s
vt=dv/vc                                                              # v-tilde
g=G*Mtot*Msun/s**2                                                    # m/s^2
with np.errstate(invalid='ignore'):
    mh=np.nanmean(np.vstack([d['mh1'],d['mh2']]),axis=0)
vtan=4.74047*np.hypot(d['pmra1'],d['pmde1'])*dist/1000.0               # km/s (thin-disk proxy)
ms=(MG1>3)&(MG1<14)&(MG2>3)&(MG2<14)
good=np.isfinite(vt)&np.isfinite(g)&np.isfinite(mh)&ms&(dv/np.maximum(edv,1e-9)>3)&(s_AU>50)
vt,g,mh,vtan=vt[good],g[good],mh[good],vtan[good]

print("="*80); print("WIDE-BINARY METALLICITY-IV TEST (Gaia DR3, d<50pc self-join)"); print("="*80)
print(f"  pairs: {len(d)} queried -> {good.sum()} after cuts (MS, [M/H], dv-SNR>3, s>50AU)")
lo=g<1e-9; hi=g>1e-8
print(f"  low-acceleration (g<1e-9): {lo.sum()}   Newtonian control (g>1e-8): {hi.sum()}")

def report(name, reg, arr_mh, arr_vt):
    if reg.sum()<30: print(f"\n  [{name}] N={reg.sum()} -- UNDERPOWERED"); return None
    mhr,vtr=arr_mh[reg],arr_vt[reg]
    print(f"\n  [{name}] N={reg.sum()}  median v~={np.median(vtr):.3f}  frac(v~>1)={np.mean(vtr>1):.3f}")
    for e0,e1 in [(-3,-0.4),(-0.4,0.0),(0.0,3)]:
        b=(mhr>e0)&(mhr<=e1)
        if b.sum()>=8: print(f"     [Fe/H]({e0:+.1f},{e1:+.1f}] N={b.sum():4d}  med v~={np.median(vtr[b]):.3f}"
                             f"  frac(v~>1)={np.mean(vtr[b]>1):.3f}")
    sl,_,r,p,se=stats.linregress(mhr,vtr)
    print(f"     slope d(v~)/d[Fe/H] = {sl:+.3f} +/- {se:.3f}  (p={p:.3f})")
    return sl,se,p

print("\n"+"-"*80); print("SANITY: control must be flat in [Fe/H] and have med v~<1, ~no v~>1 tail")
report("HIGH-g control", hi, mh, vt)
print("\n"+"-"*80); print("THE TEST (contamination predicts NEGATIVE slope; gravity predicts ZERO):")
res=report("LOW-g anomaly regime", lo, mh, vt)
# confound: thin-disk cut (metallicity correlates with age/kinematics)
thin=lo&(vtan<40)
report("LOW-g, thin-disk (v_tan<40)", thin, mh, vt)

print("\n"+"="*80); print("VERDICT"); print("="*80)
if res is None:
    print("  UNDERPOWERED at d<50pc.")
else:
    sl,se,p=res
    npred=-0.20   # crude contamination prediction: dv~ per dex ~ -0.2 (multiplicity ~doubles to [Fe/H]=-1)
    print(f"  low-g slope = {sl:+.3f}+/-{se:.3f}. Contamination prediction ~ {npred:+.2f}/dex.")
    if p<0.05 and sl<0:   print("  => NEGATIVE & significant: consistent with hidden-companion contamination.")
    elif se>0.15:         print(f"  => UNDERPOWERED to distinguish 0 from {npred}: need ~{int(((se/0.05)**2)*lo.sum())} low-g pairs for 0.05 precision.")
    else:                 print("  => slope consistent with 0 and inconsistent with contamination => favors GRAVITATIONAL origin.")
print("  NOTE: prototype-scale (d<50pc); GSP-Phot [M/H] covers ~20% of nearby stars; the crude")
print("  mass-luminosity relation and contamination transfer are order-of-magnitude. See the note.")
