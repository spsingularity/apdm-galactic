"""§9 IMMEDIATE STEP, taken literally: write V(psi)=a psi^2 + c psi^6 into an actual
scalar action, compute T_mu_nu, and check the ASSERTED w(psi)=(psi-1)/2.

Last turn (test_covariant_landau_stability.py) did the STABILITY half of the gate with
schematic fields. This turn does the other half the handout demands (§9): promote the
Landau potential to a genuine covariant action and derive w(psi) from T_mu_nu -- then
compare to the handout's asserted map w=(psi-1)/2 (equivalently psi=1+2w).

Canonical scalar (the minimal covariant completion):
    S = integral d^4x sqrt(-g) [ -1/2 (d psi)^2 - V(psi) ],   V(psi)=a psi^2 + c psi^6.
Homogeneous FRW:  rho = 1/2 psidot^2 + V,   p = 1/2 psidot^2 - V,   w = p/rho,
                  EoM:  psiddot + 3 H psidot + V'(psi) = 0.

The potential is FIXED, not chosen, by the handout's own three inputs:
  (i)   even in psi        (the D-duality psi -> -psi of §6),
  (ii)  near-tricritical   (b -> 0, so no psi^4 term),
  (iii) the two BROKEN phases are the minima at psi = +1 (matter) and psi = -1 (vacuum),
        while psi = 0 is the UNSTABLE critical point (§5.2, §6).
  => minima at psi^4 = -a/(3c) = 1  =>  a = -3c.  Take c=1:  V(psi)=psi^6 - 3 psi^2.

FINDING (previewed): the field-theory w(psi) does NOT reproduce the asserted (psi-1)/2.
A static field at EITHER minimum gives w=-1; a matter face (w=0) needs OSCILLATION about
a minimum with zero energy floor, but the even double-well forces V(+1)=V(-1), so one
cannot get w=0 at one minimum and w=-1 at the other. The D-symmetry that FORCES the even
potential (§6) is incompatible with the asymmetric two-face physics it is meant to carry.
This makes §8.2's "branch-tracking constructed, not forced" precise, and shows the gate,
taken literally as a single fundamental scalar, does not pass. Reported honestly.

Refs: M. S. Turner, PRD 28, 1243 (1983) -- oscillating scalar has <w>=(n-2)/(n+2).
"""
import numpy as np
from scipy.integrate import solve_ivp

c = 1.0
a = -3.0 * c                      # pins minima to psi = +-1 (the broken phases)
V   = lambda p: a * p**2 + c * p**6
Vp  = lambda p: 2 * a * p + 6 * c * p**5
Vpp = lambda p: 2 * a + 30 * c * p**4

# ----------------------------------------------------------------------------
# PART 1 -- the potential shape: psi=0 is the UNSTABLE critical point (from V itself)
# ----------------------------------------------------------------------------
print("=" * 76)
print("PART 1  V(psi)=psi^6-3psi^2 : is psi=0 really the unstable critical point?")
print("=" * 76)
print("   psi     V(psi)     V''(psi)    character")
for p in [-1.0, 0.0, 1.0]:
    ch = "MAX  -> UNSTABLE critical point" if Vpp(p) < 0 else "MIN  -> stable broken phase"
    print(f"   {p:+.1f}    {V(p):+.2f}      {Vpp(p):+.2f}      {ch}")
print("""
   => The actual potential realizes the handout's assertion RIGOROUSLY: psi=0 is a
      local maximum (V''<0), the two broken phases psi=+-1 are the degenerate minima.
      This UPGRADES §5.2/§6 from 'the symmetric phase is pathological' (hand-wave) to
      a theorem about V. [STRENGTHENS §5.2, §6]""")

# ----------------------------------------------------------------------------
# PART 2 -- w from T_mu_nu for a STATIC field: gives w=-1 at BOTH minima
# ----------------------------------------------------------------------------
print("=" * 76)
print("PART 2  w(psi) from T_mu_nu, static field (psidot=0):  rho=V, p=-V")
print("=" * 76)
print("   phase        psi    field w = p/rho    ASSERTED w=(psi-1)/2    match?")
for p, name in [(+1.0, "matter"), (-1.0, "vacuum")]:
    rho, pr = V(p), -V(p)
    w_field = pr / rho
    w_assert = (p - 1) / 2
    print(f"   {name:9s}   {p:+.0f}      {w_field:+.2f}              {w_assert:+.2f}"
          f"                {'YES' if abs(w_field-w_assert) < 1e-9 else 'NO'}")
print("""
   => A field frozen at a minimum is vacuum-like (w=-1) REGARDLESS of which minimum.
      The asserted map wants w=0 (matter) at psi=+1 -- NOT produced. So w=(psi-1)/2 is
      equation-of-state BOOKKEEPING (a recoding of w), not the dynamics of a field in
      V(psi). [The make-or-break identification is NOT reproduced by the naive action]""")

# ----------------------------------------------------------------------------
# PART 3 -- the matter face needs OSCILLATION, and the offset decides w
# ----------------------------------------------------------------------------
print("=" * 76)
print("PART 3  Oscillating field about psi=+1 : <w> is set by the energy FLOOR V_min")
print("=" * 76)

def avg_w(V_off, disp=0.05, t1=4000.0):
    """time-averaged w for psi oscillating about the psi=+1 minimum in a matter-dom
    background (a ~ t^{2/3}, H = 2/3t), with the potential floor shifted to V_min=V(1)-V_off."""
    Vt  = lambda p: (a * p**2 + c * p**6) - V_off
    Vtp = lambda p: 2 * a * p + 6 * c * p**5
    rhs = lambda t, y: [y[1], -3 * (2.0 / (3.0 * t)) * y[1] - Vtp(y[0])]
    sol = solve_ivp(rhs, [1.0, t1], [1.0 + disp, 0.0], rtol=1e-9, atol=1e-12,
                    dense_output=True, max_step=1.0)
    tt = np.linspace(0.5 * t1, t1, 20000)
    p, pd = sol.sol(tt)
    KE, PE = 0.5 * pd**2, Vt(p)
    return np.mean(KE - PE) / np.mean(KE + PE)

print("   V_min (floor)    <w> (oscillation)    behaves as")
for V_off in [V(1.0), V(1.0) - 0.5, V(1.0) - 2.0]:
    Vmin = V(1.0) - V_off
    w = avg_w(V_off)
    beh = "pressureless MATTER (w=0)" if abs(w) < 0.05 else "cosmological-constant-like (w=-1)"
    print(f"   {Vmin:+.2f}           {w:+.3f}               {beh}")
print("""
   => Same oscillating field: <w>=0 (matter) ONLY when the floor V_min=0; ANY positive
      floor drives <w> -> -1. The energy OFFSET at the minimum is what distinguishes
      'matter' from 'dark energy' -- not the value of psi.""")

# ----------------------------------------------------------------------------
# PART 4 -- the structural obstruction: even V(psi) cannot carry both faces
# ----------------------------------------------------------------------------
print("=" * 76)
print("PART 4  The D-symmetry obstruction")
print("=" * 76)
print(f"   even double-well forces V(+1) = V(-1) = {V(1.0):+.2f}  (the two minima are DEGENERATE)")
print("""
   To realize the two faces one needs, simultaneously,
       psi=+1 : oscillating with floor V_min = 0     -> w = 0   (matter)
       psi=-1 : frozen with floor V_min > 0          -> w = -1  (dark energy)
   but V(+1)=V(-1) makes the floors EQUAL. Set them to 0 and the vacuum face carries
   zero energy (no dark energy at all); lift them and the matter face is CC too. So an
   even (D-symmetric) potential CANNOT give a w=0 matter face AND a w=-1 DE face at once.
   The very symmetry §6 uses to FORCE the even near-tricritical form is incompatible with
   the asymmetric two-face physics it is meant to carry.

   VERDICT: The §9 gate, taken literally as ONE fundamental scalar psi with V(psi), does
   NOT pass: the asserted w(psi)=(psi-1)/2 is not the field-theory w, and the D-symmetry
   obstructs a single-field realization. What IS earned: (1) psi=0 is rigorously the
   unstable critical point (Part 1) [STRENGTHENS §5.2/§6]; (2) w=(psi-1)/2 is a STRUCTURAL
   recoding of the equation of state, not field dynamics -- so V(psi) is a Landau FREE
   ENERGY over the order-parameter coordinate, NOT a scalar-field potential. This makes
   §8.2 precise and DEMOTES §6's 'the phase structure IS the DE sector' to 'the phase
   structure is a diagnostic OF the DE sector'. The covariant carriers remain the BK
   superfluid (DM face) and the SEDE/Barrow horizon sector (DE face); psi just LABELS
   their combined w. [HONEST DOWNGRADE + two salvaged results]""")
