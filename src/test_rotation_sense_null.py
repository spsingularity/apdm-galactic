"""KILL TEST for IDEA 2 -- the rotation-sense null test (co- vs counter-rotating tracers).

Proposed smoking gun: MOND-as-a-law and collisionless CDM are rotation-blind (the RAR is identical
for co- and counter-rotating components in the same potential). A SUPERFLUID halo carries its
angular momentum as a quantized-vortex lattice and has a bulk irrotational flow v_s=grad(theta)/m;
if baryons couple to the RELATIVE velocity w = v_s - v_b, then co- and counter-rotating tracers in
the same galaxy would sit on DIFFERENT RARs -- a clean law-vs-medium discriminant in existing data
(MaNGA counter-rotators, NGC 4550-type systems).

Kill test (theory + order-of-magnitude, decided by the numbers below): does the BK EFT actually
generate an OBSERVABLE relative-velocity force?  Three sub-claims:

  A. LEADING ORDER is exactly rotation-blind. BK's coupling L = -alpha (Lambda/M_Pl) theta rho_b
     couples the phonon to baryon DENSITY, not current => F = -alpha(Lambda/M_Pl) grad(theta),
     independent of v_b. Co and counter feel the identical force. [the naive smoking gun is absent]
  B. The differential co/counter signal is CAPPED BY THE HALO SPIN. The vortex lattice is ultra-fine
     => the superfluid coarse-grains to solid-body flow at v_s ~ lambda*v_c (spin lambda~0.035), not
     v_c. So |w_co| and |w_counter| differ only by ~2*lambda, whatever the coupling.
  C. Both channels that could break the null are dead: the conservative relative-velocity operator
     rho_b*w^2 is symmetry-allowed but absent in minimal BK and would itself broaden the RAR if O(1);
     the dissipative vortex mutual friction is ~ (rho_n/rho) ~ (T/T_c)^4, negligible.

PASS bar (pre-registered): a clean smoking gun needs a differential co/counter RAR offset that is
(i) present at leading order OR (ii) >~ the RAR scatter (0.11 dex) at a coupling the theory actually
contains. Verdict decided by the numbers.
"""
import numpy as np

# ---- constants (SI for the vortex lattice) ----
hbar = 1.054571817e-34      # J s
c    = 2.99792458e8         # m/s
eV   = 1.602176634e-19      # J
kpc  = 3.0856775814913673e19  # m
AU   = 1.495978707e11       # m

# ---- repo parameter point (paper sec.3) ----
m_eV   = 1.0                # DM boson mass [eV]
m_kg   = m_eV*eV/c**2       # [kg]
v_c    = 150e3             # HSB rotation speed [m/s]
T_over_Tc = 2e-4           # paper sec.3 (galaxy: T_DM/T_c)
lam    = 0.035             # halo spin parameter (Bullock+2001 median)
R_halo = 30*kpc            # halo scale radius [m]

print("="*80)
print("IDEA 2 KILL TEST: co- vs counter-rotating RAR offset in BK superfluid DM")
print("="*80)
print(f"  point: m={m_eV} eV, v_c={v_c/1e3:.0f} km/s, T/T_c={T_over_Tc:.0e}, spin lambda={lam}, "
      f"R={R_halo/kpc:.0f} kpc")

print("\n" + "-"*80)
print("A. LEADING-ORDER FORCE IS ROTATION-BLIND (the naive smoking gun does not exist)")
print("-"*80)
print("  BK coupling  L_int = -alpha (Lambda/M_Pl) theta * rho_b   couples phonon to baryon DENSITY.")
print("  => F/m_b = -alpha (Lambda/M_Pl) grad(theta) = function of the LOCAL superfluid field only.")
print("  F(v_b) - F(-v_b) = 0  identically: a co- and a counter-rotating baryon at the same point")
print("  feel the SAME MOND force.  So the leading RAR is identical for both -> NO LO signal.")
print("  (This is robust to the exact coupling: any coupling to rho_b, not the current j_b=rho_b v_b,")
print("   gives a velocity-independent force. The current-coupling needed for a signal is sub-leading.)")

print("\n" + "-"*80)
print("B. THE VORTEX LATTICE IS ULTRA-FINE -> superfluid coarse-grains to v_s ~ lambda*v_c")
print("-"*80)
kappa = 2*np.pi*hbar/m_kg               # circulation quantum h/m [m^2/s]
v_s   = lam*v_c                         # bulk superfluid flow (spin-limited) [m/s]
Omega = v_s/R_halo                      # angular velocity [1/s]
n_v   = 2*Omega/kappa                   # vortex areal density [1/m^2]
ell_v = 1.0/np.sqrt(n_v)               # inter-vortex spacing [m]
N_v   = n_v*np.pi*R_halo**2            # total number of vortices
print(f"  circulation quantum kappa = h/m      = {kappa:.1f} m^2/s   (huge, since m is tiny)")
print(f"  bulk flow  v_s ~ lambda*v_c          = {v_s/1e3:.1f} km/s   (<< v_c={v_c/1e3:.0f} km/s)")
print(f"  vortex spacing  ell_v = 1/sqrt(n_v)  = {ell_v:.2e} m = {ell_v/AU:.3f} AU  (<< kpc)")
print(f"  total vortices  N_v                  = {N_v:.1e}")
print(f"  => spacing {ell_v/kpc:.1e} kpc << galactic: the lattice coarse-grains to SMOOTH solid-body")
print(f"     rotation at v_s ~ {v_s/1e3:.0f} km/s. Per-vortex Magnus kicks average out.")
diff_ceiling = 2*lam                     # |w_counter| - |w_co| over v_c
print(f"  differential ceiling: |w_counter|-|w_co| = 2*v_s = 2*lambda*v_c = {diff_ceiling:.3f}*v_c")
print(f"     -> co and counter differ by only ~{100*diff_ceiling:.0f}% in relative speed, NOT 200%.")

print("\n" + "-"*80)
print("C. BOTH NULL-BREAKING CHANNELS ARE DEAD")
print("-"*80)
# (C1) conservative relative-velocity operator rho_b * w^2 : symmetry-allowed, absent in minimal BK.
# Even at MAXIMAL O(1) coupling (force fully modulated by w), the co/counter differential is bounded
# by how much |w| changes: d a / a ~ n * (2 v_s / w) with w~v_c  ->  ~ 2*lambda per power of w.
eps_max = 2*lam                          # fractional force differential at O(1), linear-in-w coupling
dex_max = np.log10(1+eps_max)            # as a coherent RAR-normalization offset
RAR_scatter = 0.11                       # dex (Lelli+2016 / McGaugh+2016)
print(f"  C1 conservative rho_b*w^2 operator (Galilean-allowed, NOT in minimal BK):")
print(f"     max differential force (O(1) coupling) ~ 2*lambda = {eps_max:.3f} = {dex_max:.3f} dex")
print(f"     vs RAR scatter {RAR_scatter} dex -> even the CEILING is below the per-point scatter, and")
print(f"     it requires a coefficient the theory does not contain. An O(1) w-coupling would ALSO make")
print(f"     the RAR depend on v_c (30-300 km/s across the sample) -> broaden it beyond the observed")
print(f"     {RAR_scatter} dex. The RAR's own tightness bounds this coupling to be small: self-defeating.")
rho_n_over_rho = T_over_Tc**4            # phonon normal-fluid fraction (~T^4 at T<<T_c)
print(f"  C2 dissipative vortex mutual friction ~ rho_n/rho ~ (T/T_c)^4 = {rho_n_over_rho:.1e}")
print(f"     -> the normal (dissipative) component is ~{-np.log10(rho_n_over_rho):.0f} orders down;")
print(f"        mutual-friction drag is dynamically negligible over any number of orbits.")

print("\n" + "="*80)
print("VERDICT")
print("="*80)
print("  A: leading-order force is EXACTLY rotation-blind -> the clean smoking gun does NOT exist.")
print(f"  B: superfluid rotates at ~{v_s/1e3:.0f} km/s (spin-limited), not v_c -> co/counter differ by ~{100*diff_ceiling:.0f}%,")
print("     not the naive 200%; the '2 v_c contrast' picture is wrong.")
print(f"  C: the only null-breaking operators are (C1) absent-and-self-limiting (ceiling {dex_max:.3f} dex < "
      f"{RAR_scatter} dex,")
print(f"     and a large version spoils the RAR) and (C2) dissipatively dead ((T/T_c)^4~{rho_n_over_rho:.0e}).")
print("  => IDEA 2 is NOT a clean smoking gun. The rotation-sense null is a real, clean PREDICTION of")
print("     BK (a genuine consistency point: BK must show NO co/counter offset), but it does NOT")
print("     discriminate BK from MOND/CDM, because all three predict the same (null) at observable")
print("     precision. [REAL][HIGH] for the LO null; [MODEL][MED] for the suppressed residual.")
print("  Salvage: the residual differential, if ever seen at >0.03 dex, would falsify BK+MOND+CDM at")
print("     once (none predicts it) -> keep as a cheap MaNGA null-test to run, not a positive prediction.")
