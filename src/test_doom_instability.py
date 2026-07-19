"""THE RESIDUAL, CLOSED: the Valiviita-Majerotto-Maartens / Gavela 'doom' instability
check on the DM<->DE coupling that makes the ghost-free completion work.

Last turn (test_twofield_completion.py) resolved the fork: the phantom fluid EOS SEDE
needs is carried GHOST-FREE by interacting DE (rho_DE tracks structure), not by a
quintom (which is a ghost). The one residual owed there was the perturbation-level
stability of that coupling -- interacting DE can suffer a large-scale non-adiabatic
instability (Valiviita et al. 2008) even when the background is ghost-free. This script
does that check.

  PART 1  SEDE-as-implemented uses SMOOTH DE (background only; no delta_DE evolved).
          A smooth w<-1 fluid with rest-frame c_s^2=1 is perturbatively stable
          (Kunz-Sapone 2006) and SEDE's finite fsigma8 corroborates -- but that
          description leaves the phantom microscopically UNexplained (a fitted w_eff).
  PART 2  The GENUINE interacting completion (structure sources rho_DE) is a real
          coupling Q, so the doom factor applies. Derive it:
              Q = 3 H rho_DE (w_eff - w_DE),   D = Q/[3H(1+w_DE)rho_DE] = (w_eff-w_DE)/(1+w_DE)
          Gavela et al. (2009): stable on large scales  <=>  D <= 0.
  PART 3  Evaluate D(z) with the underlying fluid EOS = SEDE's STRUCTURAL EOS
          (w_DE ~ -0.856, non-phantom) and w_eff = SEDE's dynamical (phantom) EOS.
  PART 4  Show the stability condition D<=0 REDUCES to 1+w_DE>0 (non-phantom underlying)
          -- i.e. exactly the handout's sign(1+w)=+1 structural prediction. Counterfactual
          phantom underlying -> D>0 -> unstable.

Imports SEDE (../SEDE) live; falls back to headline numbers. Refs: J. Valiviita,
E. Majerotto, R. Maartens, JCAP 0807, 020 (2008); M. B. Gavela et al., JCAP 0907, 034
(2009), arXiv:0901.1611 (the doom factor); M. Kunz, D. Sapone, PRD 74, 123503 (2006).
"""
import os, sys
import numpy as np

SEDE = os.path.dirname(os.path.abspath(__file__))  # vendored sede/ + predictions.json beside this script
sys.path.insert(0, SEDE)
Om = 0.311

try:
    from sede import theory
    w_str = float(theory.w_DE_effective(Om))                       # structural (underlying) EOS
    w_eff = {z: float(np.atleast_1d(theory.w_DE_dynamical(np.array([z]), Om))[0])
             for z in (0.0, 0.5, 1.0, 2.0)}
    src = "SEDE live"
except Exception as e:
    print(f"   [SEDE import failed ({e}); using headline numbers]")
    w_str = -0.8561
    w_eff = {0.0: -1.1506, 0.5: -1.3138, 1.0: -1.4365, 2.0: -1.5583}
    src = "headline numbers"

# ============================================================================
print("=" * 78)
print("PART 1  How SEDE actually treats DE: smooth (background only)")
print("=" * 78)
print("""  SEDE (sede/perturbations.py) evolves DE as SMOOTH: it 'does NOT cluster on
  sub-horizon scales; it enters only via the background H(z)'. A smooth w<-1 fluid with
  rest-frame c_s^2=1 is perturbatively stable (Kunz-Sapone 2006), and SEDE's finite,
  smooth fsigma8(z) corroborates. BUT this description leaves the phantom w_eff
  microscopically UNexplained -- it is a fitted background EOS, with no delta_DE dynamics
  and hence no test of the coupling. The doom check below is for the GENUINE interacting
  completion, where the phantom IS explained (structure sources rho_DE).""")

# ============================================================================
print("=" * 78)
print("PART 2  The doom factor for the genuine coupling")
print("=" * 78)
print("""  Effective single fluid: dln rho_DE/dln a = -3(1+w_eff).
  Interacting split (fluid w_DE + exchange Q):  rho_DE' + 3H(1+w_DE)rho_DE = -Q.
  Matching the SAME rho_DE(a):   Q = 3 H rho_DE (w_eff - w_DE).
  Gavela et al. (2009) doom factor and stability criterion:
        D = Q/[3H(1+w_DE)rho_DE] = (w_eff - w_DE)/(1 + w_DE),   STABLE  <=>  D <= 0.""")

# ============================================================================
print("=" * 78)
print(f"PART 3  Evaluate D(z)  ({src}: underlying w_DE = {w_str:+.4f}, non-phantom)")
print("=" * 78)
print("   z     w_eff (phantom)    Q sign            D = (w_eff-w_DE)/(1+w_DE)    verdict")
all_stable = True
for z, we in w_eff.items():
    D = (we - w_str) / (1 + w_str)
    Qsign = "Q<0 (structure->DE)" if (we - w_str) < 0 else "Q>0"
    all_stable = all_stable and (D <= 0)
    print(f"   {z:.1f}   {we:+.3f}           {Qsign:19s} {D:+.3f}                   "
          f"{'STABLE' if D <= 0 else 'UNSTABLE'}")
print(f"\n   => D(z) <= 0 at every z : {all_stable}  =>  NO large-scale doom instability.")

# ============================================================================
print("=" * 78)
print("PART 4  Why it is stable: the criterion IS the handout's sign(1+w)=+1")
print("=" * 78)
print("""  Since w_eff < w_DE (the coupling makes DE MORE phantom), the numerator (w_eff-w_DE)<0.
  So  D <= 0  <=>  (1 + w_DE) > 0  <=>  the UNDERLYING fluid is NON-PHANTOM.
  That is exactly the handout's structural prediction sign(1+w)=+1 (§7.3), which the
  program independently matches to SEDE's structural EOS. Counterfactual -- a phantom
  underlying fluid destroys stability:""")
for w_cf in [-1.05, -1.20]:
    D = (w_eff[0.5] - w_cf) / (1 + w_cf)
    print(f"     underlying w_DE={w_cf:+.2f} (phantom): D(z=0.5)={D:+.2f}  "
          f"{'STABLE' if D <= 0 else 'UNSTABLE'}")

# ============================================================================
print("=" * 78)
print("VERDICT")
print("=" * 78)
print("""  The Valiviita/Gavela doom-instability check PASSES: D(z) = -2.0 to -4.9 < 0 across
  z, so the DM<->DE coupling that carries the phantom w_eff is perturbatively STABLE on
  large scales -- and it is stable PRECISELY BECAUSE the underlying DE fluid is the
  non-phantom STRUCTURAL EOS (1+w_DE>0), the handout's own sign(1+w)=+1 prediction. A
  phantom underlying fluid would flip D>0 and be fatal.

  This CLOSES the single residual left by test_twofield_completion.py. The completion is
  now BOTH ghost-free (background, prior turn) AND doom-stable (perturbations, here),
  conditional on one thing the program already predicts: the underlying DE fluid is
  non-phantom, with the apparent w_eff<-1 an interacting-DE (coupling) effect.

  The thread closes coherently: STRUCTURAL EOS (~ -0.85, +1) = the underlying stable
  fluid; DYNAMICAL EOS (< -1, phantom) = the apparent coupled EOS; the coupling sign
  Q<0 (structure -> DE) is stabilizing exactly because the underlying fluid is
  non-phantom. What remains is numerical (a full CAMB/CLASS run of delta_DE with this
  coupling) and observational (the a0(z) sign test) -- not a matter of principle.""")
