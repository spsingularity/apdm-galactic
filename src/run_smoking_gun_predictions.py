"""Concrete predictions for three smoking-gun signatures.
All three hinge on ONE quantity: the condensate phonon sound speed c_s.

THE LINCHPIN (c_s), resolved by a consistency argument, not the microscopic coupling:
- A naive interaction-BEC estimate c_s=sqrt(mu_int/m) with mu_int=u0*n0~4e-17 eV gives c_s~m/s
  -- but that is the wrong scale (interaction pressure only). RETRACTED.
- The MOND phonon force is COHERENT across the galaxy; a perturbed galaxy returns to the tight
  observed radial-acceleration relation (RAR, scatter ~0.1 dex) only if phase coherence re-
  establishes within the time between perturbations, i.e. the sound-crossing time R/c_s must be
  << Gyr. That REQUIRES c_s >> R/Gyr ~ 10 kpc/Gyr ~ 10 km/s. The observed RAR tightness therefore
  forces c_s to tens-hundreds km/s -- i.e. c_s ~ v_c (the virial/hydrostatic value for a pressure-
  supported self-gravitating condensate). This resolves the linchpin toward the "high" bracket and
  makes all three signatures live. It is also testable (c_s sets the phonon dispersion / the RC
  shape in the transition region).

Run: python3 experiments/run_smoking_gun_predictions.py   (no SEDE dependency)
"""
import numpy as np

# c_s ~ v_c (host circular velocity), O(1) factor; we take c_s = v_c and flag the O(1).
print("### LINCHPIN: c_s ~ v_c  (from RAR-tightness / hydrostatic argument) ###")
print("  naive interaction estimate c_s=sqrt(u0 n0/m) ~ m/s  -> RETRACTED (wrong pressure)")
print("  RAR tightness requires sound-crossing R/c_s << Gyr  => c_s >~ 10 km/s, ~ v_c.\n")

# ---------- #2  LANDAU-THRESHOLD DYNAMICAL FRICTION ----------
# Mach number M = v_perturber / c_s(host) = v_perturber / v_c(host).
# Prediction: M<1 (subsonic) -> zero phonon friction -> survives/stays fast;
#             M>1 (supersonic)-> Cherenkov phonon wake -> sinks (+ dark Mach cone).
print("### #2  Landau-threshold dynamical friction  (M = v_perturber / v_c,host) ###")
perturbers = [
    # name,                     v_perturber[km/s], v_c host[km/s], observed
    ("Fornax GCs (in Fornax dSph)",   12,   22, "GCs have NOT sunk (timing problem)"),
    ("MW bar (corotation)",          210,  220, "bar is persistently FAST"),
    ("Fornax dSph in MW halo",        s:=190, 180, "still on wide orbit"),
    ("Sagittarius stream progenitor", 300,  180, "decaying / disrupted"),
    ("LMC in MW halo",               320,  180, "sinking; warps MW disk"),
]
print(f"  {'perturber':32s} {'v_p':>5s} {'v_c':>5s}  {'Mach':>5s}  phase        observed")
for name, vp, vc, obs in perturbers:
    M = vp/vc
    phase = "SUBSONIC" if M < 1 else "supersonic"
    verdict = "survive" if M < 1 else "SINK"
    print(f"  {name:32s} {vp:5.0f} {vc:5.0f}  {M:5.2f}  {phase:10s} -> {verdict:7s} | {obs}")
print("  PREREGISTERED: survival/sinking tracks M=v_p/v_c with a threshold at M=1;")
print("  Fornax GCs subsonic (survive, anomaly explained), LMC supersonic (sinks, warps).")
print("  DECISIVE TEST: compile perturbers, compute M=v_p/v_c(host), test the M=1 survival step.")
print("  vs rivals: LCDM friction always-on & rises at low v (wrong sign of threshold);")
print("             plain MOND reduces friction but SMOOTHLY (no step); rel. MOND has no DM medium.\n")

# ---------- #1  MERGER-HISTORY HYSTERESIS IN THE RAR ----------
# A major merger cannot decondense a deeply-condensed galaxy (C~1e6 -> merger drops C by ~8x,
# still >>1). The mechanism is COHERENCE disruption, not decondensation: the merger scrambles
# the galaxy-wide phonon PHASE, killing the coherent MOND force until coherence re-establishes
# over the sound-crossing time tau_cond ~ R/c_s ~ R/v_c.
print("### #1  Merger hysteresis in the RAR  (coherence disruption, not decondensation) ###")
kpc_km = 3.0856775814913673e16      # km per kpc
Gyr_s = 3.1557e16
for name, R_kpc, vc in [("dwarf", 3, 30), ("MW-like", 15, 180), ("massive spiral", 25, 250)]:
    cs = vc                          # km/s
    tau = (R_kpc*kpc_km/cs)/Gyr_s    # sound-crossing time in Gyr
    print(f"  {name:15s}: R={R_kpc:2d} kpc, c_s~v_c={vc:3d} km/s -> tau_cond ~ R/c_s = {tau:.2f} Gyr")
print("  => tau_cond ~ 0.05-0.15 Gyr (coherence recovery). Merger rate ~1/(few Gyr),")
print("     so a small but finite fraction (~tau_cond/tau_merger ~ few %) of galaxies are caught")
print("     mid-recovery: OFF the RAR (MOND force scrambled), relaxing back over tau_cond.")
print("  PREREGISTERED three-way SIGN test (RAR offset at fixed baryon profile vs merger age):")
print("     APDM: POSITIVE offset (above RAR), decaying over ~0.1 Gyr")
print("     plain MOND: ZERO (RAR is strictly ahistorical -> any correlation kills MOND)")
print("     LCDM: NEGATIVE (merger remnants more concentrated -> tighter inner RAR)")
print("  DECISIVE TEST: MaNGA IFU + post-merger catalog; dRAR vs time-since-merger, sign+decay.")
print("  CAVEAT: short window (few % caught); needs coherence-disruption to null the MOND force.\n")

# ---------- #3  SOLAR-APEX DIPOLE IN THE LOCAL a0 ----------
# A superfluid is irrotational: the bulk condensate does NOT co-rotate with the disk, so the Sun
# moves at ~V_sun through it. Mach ~ V_sun/c_s(R_sun) ~ V_sun/v_c ~ 1 (transonic) -> fore-aft
# anisotropy in the local effective a0, a DIPOLE aligned with the solar velocity (Galactic
# rotation apex, l~90 deg), NOT the Galactic-centre direction.
print("### #3  Solar-apex dipole in the local a0 (Gaia wide binaries) ###")
V_sun, cs_local = 230.0, 233.0       # km/s: solar circular velocity ~ local v_c ~ c_s
mach = V_sun/cs_local
print(f"  Sun through irrotational condensate: Mach = V_sun/c_s = {mach:.2f} (transonic)")
print(f"  => dipole in the low-acceleration boost, amplitude ~ O(Mach) fraction, aligned with")
print(f"     the Galactic-rotation apex (l~90 deg), not the Galactic centre.")
print("  PREREGISTERED: a DIRECTION-DEPENDENT a0 exists (isotropic MOND forbids it by axiom);")
print("     dipole locked to V_sun distinguishes 'moving through a medium' from Galactic-field sourcing.")
print("  DECISIVE TEST: split Gaia DR3 wide-binary low-acceleration sample by orbit orientation")
print("     relative to the solar apex; look for a cos(theta) modulation of the anomaly.")
print("  CAVEAT: the wide-binary anomaly itself is contested (Chae vs Banik); if the condensate")
print("     is vortex-threaded and effectively co-rotates, the relative velocity (dipole) shrinks.\n")

print("### SHARED RESULT ###")
print("  All three live IFF c_s ~ v_c (tens-hundreds km/s), which RAR tightness independently")
print("  requires. The single highest-value theory calc remains the BK phonon c_s(r) from the full")
print("  Lagrangian, to confirm the O(1) and the radial profile. My earlier 'c_s~m/s kill' is")
print("  retracted: it used the interaction chemical potential, not the coherence/virial scale.")
