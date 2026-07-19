"""§9 DOWNSTREAM tasks -- the three original numerical tasks, now framed as
computations on V(psi)/the phase structure:

  TASK A  sign(1+w_SEDE):  the handout (§7.3) predicts +1 (non-phantom) structurally.
          Verify against the EXISTING SEDE derivations (../SEDE, imported live).
  TASK B  SIDM-preferred point in the near-tricritical window:  does the BK/APDM SIDM
          benchmark map to the physical psi-range near psi->+1 with the observed sigma(v)?
  TASK C  Thermalization clock vs the critical surface psi=0:  is the condensation time
          of psi (the phase transition) the SAME object as the cosmological psi=0 surface
          (rho_DM = rho_DE), i.e. does the thermalization clock track structure formation?

These build directly on the two action experiments (test_covariant_landau_stability.py,
test_landau_action_wpsi.py). Task A is quantitative (SEDE numbers); B and C are honest
mappings that say what is and is not established.
"""
import os, sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # vendored sede/ beside this script

# ============================================================================
# TASK A -- sign(1+w) : handout §7.3 (+1) vs the actual SEDE derivations
# ============================================================================
print("=" * 76)
print("TASK A  sign(1+w_SEDE):  handout §7.3 says +1 (non-phantom). Check SEDE.")
print("=" * 76)
Om = 0.311
try:
    from sede import theory
    w_struct = theory.w_DE_effective(Om)          # structural (order-parameter) EOS
    w_alg    = theory.w_DE_algebraic(Om)          # structural, algebraic approx
    w_fluid0 = theory.w_DE_fluid(Om)              # dynamical fluid EOS at z=0
    w_dyn = {z: float(np.atleast_1d(theory.w_DE_dynamical(np.array([z]), Om))[0])
             for z in (0.0, 0.5, 1.0, 2.0)}
    have_sede = True
except Exception as e:  # pragma: no cover -- fall back to SEDE/predictions.json values
    print(f"   [SEDE import failed: {e}; using SEDE/predictions.json headline numbers]")
    w_struct, w_alg, w_fluid0 = -0.8561, -0.8495, -1.1506
    w_dyn = {0.0: -1.1506, 0.5: -1.3138, 1.0: -1.4365, 2.0: -1.5583}
    have_sede = False

print("\n  (A1) SEDE STRUCTURAL EOS (derivative w.r.t. the order parameter x=D^2):")
for name, w in [("w_DE_effective", w_struct), ("w_DE_algebraic", w_alg)]:
    print(f"        {name:16s} = {w:+.4f}   1+w = {1+w:+.3f}   sign(1+w) = {'+1' if 1+w>=0 else '-1'}")
print("        => matches handout §7.3: sign(1+w) = +1  (NON-PHANTOM).")

print("\n  (A2) SEDE DYNAMICAL FLUID EOS (the object that governs H(z); d/d ln a):")
print(f"        w_DE_fluid(z=0)  = {w_fluid0:+.4f}   1+w = {1+w_fluid0:+.3f}   PHANTOM")
for z, w in w_dyn.items():
    print(f"        w_dynamical(z={z:.1f}) = {w:+.4f}   1+w = {1+w:+.3f}   "
          f"{'PHANTOM  sign(1+w)=-1' if w < -1 else 'non-phantom'}")
print("""
   VERDICT (A): SEDE carries BOTH signs. The handout's sign(1+w)=+1 matches SEDE's
   STRUCTURAL EOS (~ -0.85) but is CONTRADICTED by SEDE's DYNAMICAL FLUID EOS (~ -1.15
   to -1.56, phantom at all z) -- the one that actually governs the expansion and the
   FALLING a0(z) the program uses. This is exactly the split found independently in
   test_landau_action_wpsi.py: psi=1+2w is a STRUCTURAL (order-parameter) coordinate,
   so the handout's prediction lives at the structural/bookkeeping level and matches
   SEDE's structural EOS -- NOT the gravitating fluid EOS. And it is the same phantom
   crossing the Vikman no-go (test_covariant_landau_stability.py, Part C) says a single
   ghost-free scalar cannot produce. Three routes, one conclusion. [PARTIAL: +1 holds
   structurally, FAILS dynamically]""")

# ============================================================================
# TASK B -- SIDM-preferred point vs the psi-range / near-tricritical window
# ============================================================================
print("=" * 76)
print("TASK B  Does the SIDM-preferred point sit near psi->+1 (matter phase)?")
print("=" * 76)
# The BK/APDM SIDM benchmark (README; compute_solidification.py): the condensate lives
# and self-interacts in galaxies/groups -- i.e. in the BROKEN MATTER phase, which the
# handout labels psi=+1 (w=0). Thermalization split (BK): dwarf ~1e6, group ~11, cluster ~0.013.
split = {"dwarf (~30 km/s)": 1e6, "group (~260 km/s)": 11.0, "cluster (~1500 km/s)": 0.013}
print("  BK thermalization ratio t_H/t_relax (condensed if >1) at the benchmark point:")
for env, r in split.items():
    print(f"    {env:22s} {r:>10.3g}   {'CONDENSED (psi=+1 matter phase)' if r > 1 else 'not condensed'}")
print("""
   VERDICT (B): The SIDM cross-section is fixed where the condensate EXISTS -- galaxies
   and groups -- which is precisely the broken MATTER phase psi=+1 (w=0). This is where
   last turn's STABILITY analysis also places the local field (mu>0 only as psi->+1). So
   the SIDM point maps consistently to psi->+1 with the observed sigma(v). BUT 'near-
   tricritical' is a statement about the POTENTIAL's shape (b->0), and the SIDM point --
   sitting DEEP in the broken phase, far from the critical point psi=0 -- does NOT test
   b->0. Honest reading: CONSISTENT with the phase assignment, NOT a confirmation of
   near-tricriticality. [CONSISTENT, not confirmatory]""")

# ============================================================================
# TASK C -- thermalization clock vs the cosmological critical surface psi=0
# ============================================================================
print("=" * 76)
print("TASK C  Is the psi-condensation clock the SAME object as the psi=0 surface?")
print("=" * 76)
Om, ODE = 0.315, 0.685
# cosmological psi=0  <=>  w_tot=-1/2  <=>  rho_DM=rho_DE  <=>  Om(1+z)^3 = ODE
z_eq = (ODE / Om) ** (1.0 / 3.0) - 1.0
print(f"  cosmological critical surface psi=0 (rho_DM=rho_DE):  z_eq = {z_eq:.2f}")
print(f"  BK condensation / galaxy assembly (structure formation): z ~ 1 - 3")
print(f"  ratio in (1+z):  critical (1+z)={1+z_eq:.2f}  vs  formation (1+z)~2-4")
print("""
   VERDICT (C): The two clocks are NUMERICALLY ADJACENT (order unity in 1+z) but NOT
   shown identical: rho_DM=rho_DE at z~0.3, while the condensate/structure forms at
   z~1-3. Their proximity IS the 'why now' coincidence itself (handout §1, still listed
   unexplained in §8) -- restating it, not deriving it. Claiming the thermalization
   clock and the critical surface psi=0 are 'the same object' is NOT established here;
   they coincide only to order unity, which is the coincidence to be explained, not an
   explanation. [NOT ESTABLISHED -- remains the open 'why now']""")

print("=" * 76)
print("SUMMARY of downstream tasks")
print("=" * 76)
print("""  A  sign(1+w): +1 holds for SEDE's STRUCTURAL EOS (matches §7.3), FAILS for the
       dynamical fluid EOS (phantom) -- psi=1+2w is a structural coordinate. [PARTIAL]
  B  SIDM point maps to psi=+1 (matter phase), consistent with stability; does NOT
       test near-tricriticality (b->0). [CONSISTENT, not confirmatory]
  C  psi-condensation clock vs psi=0 surface: adjacent to order unity = the 'why now'
       coincidence, not a derived identity. [NOT ESTABLISHED]""")
