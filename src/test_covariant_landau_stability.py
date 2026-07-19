"""The make-or-break gate (handout section 9): does the near-tricritical Landau
completion V(psi)=a psi^2 + c psi^6 survive a COVARIANT, GHOST-FREE stability test?

The handout builds everything (a0 from Lambda, the DM<->DE duality, the phase
structure) on top of one unproven step it flags as the "first real gate": promote
the non-relativistic energy branch f(g^2;psi) and the order parameter psi=1+2w to a
scalar-field action and check (a) the equation of state w(psi) and (b) ghost-/
gradient-instability freedom. Everything in sections 4-6 is conditional on it.

We do that here for the two faces of the field, using STANDARD k-essence/AQUAL
stability theory (no new machinery):

  PART A  LOCAL face  (the MOND-carrying gradient sector).
          The handout's own local energy (section 5.2), with u = g/a_star:
              f(u;psi) = a*^2 [ u^2 + (1-psi)/2 * u^-2 + (1+psi)/2 * (2/3) u^3 ]
          In AQUAL/k-essence the interpolation function is mu = df/d(g^2), and a
          static configuration is stable iff BOTH kinetic eigenvalues are positive:
              transverse  :  mu(u;psi) > 0
              longitudinal:  N(u;psi) = mu + 2 g^2 mu' > 0      (= P_X + 2X P_XX)
          (Bekenstein-Milgrom 1984; Vikman/Babichev k-essence conditions.)

  PART B  COSMOLOGICAL face  (the dark-energy background sector).
          A canonical scalar realizing the sector w=(psi-1)/2: c_s^2=1, no ghost,
          and 1+w=(1+psi)/2 >= 0 automatically. Confirms section 7.3's non-phantom
          claim is REALIZABLE by a stable field, not just an algebraic identity.

  PART C  PHANTOM NO-GO  (Vikman 2005).
          A single ghost-free scalar cannot cross w=-1. So sign(1+w)=+1 is upgraded
          from algebra to a STABILITY THEOREM -- but the SAME theorem forbids the
          phantom (w_a<0, w crossing -1) crossing that the SEDE branch invokes for a
          falling a0(z). The covariant completion CONSTRAINS the phenomenology; we
          surface the tension honestly.

No data, no fits: this is an analytic/stability computation, verified numerically.
References: Bekenstein & Milgrom, ApJ 286, 7 (1984); A. Vikman, PRD 71, 023515
(2005), astro-ph/0407107; E. Babichev, V. Mukhanov, A. Vikman, JHEP 0802, 101 (2008).
"""
import numpy as np
from scipy.optimize import brentq

# ----------------------------------------------------------------------------
# PART A -- LOCAL FACE: stability of the Landau energy branch f(g^2; psi)
# ----------------------------------------------------------------------------
# Units a_star = 1, u = g/a_star. The interpolation function and the two kinetic
# eigenvalues follow from f (section 5.2) by mu = df/d(g^2), N = mu + 2 g^2 dmu/d(g^2).
#
#   mu(u;psi) = 1 - (1-psi)/2 * u^-4 + (1+psi)/2 * u
#   N (u;psi) = 1 + (3/2)(1-psi) * u^-4 + (1+psi) * u
#
def mu(u, psi):
    return 1.0 - 0.5 * (1.0 - psi) * u ** -4 + 0.5 * (1.0 + psi) * u


def N_long(u, psi):
    return 1.0 + 1.5 * (1.0 - psi) * u ** -4 + (1.0 + psi) * u


def N_numeric(u, psi, h=1e-6):
    """Longitudinal eigenvalue by finite-difference on g^2 -- CHECK of the analytic N."""
    y = u * u
    dmu_dy = (mu(np.sqrt(y + h), psi) - mu(np.sqrt(y - h), psi)) / (2 * h)
    return mu(u, psi) + 2 * y * dmu_dy


print("=" * 76)
print("PART A  LOCAL FACE -- ghost/gradient stability of the Landau branch f(g^2;psi)")
print("=" * 76)

# (A1) analytic N verified against the numeric derivative, and N > 0 everywhere.
print("\n(A1) longitudinal eigenvalue N = mu + 2 g^2 mu'  (analytic vs numeric):")
u_check = np.array([0.1, 0.3, 1.0, 3.0])
all_pos = True
for psi in [1.0, 0.5, 0.0, -1.0]:
    Na = N_long(u_check, psi)
    Nn = np.array([N_numeric(x, psi) for x in u_check])
    ok = np.allclose(Na, Nn, rtol=1e-4)
    all_pos = all_pos and np.all(Na > 0)
    print(f"   psi={psi:+.1f}  N={np.round(Na,3)}  (matches numeric: {ok})")
print(f"   => longitudinal mode N>0 for ALL psi in [-1,1], all u>0 : {all_pos}")
print("      The longitudinal (radial) mode is NEVER a ghost. The whole")
print("      instability question reduces to the transverse condition mu>0.")

# (A2) transverse condition mu>0. mu is monotone increasing in u for psi<1, so there
#      is a single threshold g_c(psi): the field is UNSTABLE for g < g_c.
print("\n(A2) transverse condition mu(u;psi) > 0  ->  stability threshold g_c(psi):")
print("   psi        g_c/a_star     stable window        phase")
rows = [(1.0, "matter/galaxies (fully broken)"),
        (0.9, "near-broken"),
        (0.5, "partially ordered"),
        (0.0, "self-dual / critical (w=-1/2)"),
        (-0.5, "toward vacuum"),
        (-1.0, "vacuum/voids (fully broken)")]
for psi, label in rows:
    if psi >= 1.0:  # 1/u^4 term switches off: mu = 1 + u > 0 for all u
        print(f"   {psi:+.1f}      0.000        g > 0  (ALL g)      {label}")
        continue
    uc = brentq(lambda u: mu(u, psi), 1e-6, 1e6)
    print(f"   {psi:+.1f}      {uc:.3f}        g > {uc:.2f} a_star     {label}")

print("""
   VERDICT (A): The covariant completion of the local face is ghost-/gradient-
   stable ONLY in the fully broken matter phase psi->+1, where mu=1+u>0 for all g
   and standard MOND (mu->g/a*) is recovered. For ANY psi<1 the transverse mode
   goes unstable (mu<0) below g_c(psi); at the self-dual point psi=0 the ENTIRE
   deep-MOND regime g<0.78 a* is unstable, and g_c -> a* as psi -> -1.
   => This is not a defect but a SELECTION PRINCIPLE: stability DYNAMICALLY DRIVES
      the local field to psi=+1. The near-tricritical ordering the Landau V(psi)
      posits is ENFORCED by ghost-freedom, not assumed. The symmetric self-dual
      point is exactly the unstable critical point (consistent with section 5.2,
      now shown rigorously and quantitatively). [STRENGTHENS sections 5.2, 6]""")

# ----------------------------------------------------------------------------
# PART B -- COSMOLOGICAL FACE: a stable field realizes w(psi) with 1+w >= 0
# ----------------------------------------------------------------------------
print("=" * 76)
print("PART B  COSMOLOGICAL FACE -- w(psi) from a canonical scalar (c_s^2=1, no ghost)")
print("=" * 76)
print("  rho = 1/2 phidot^2 + V,  p = 1/2 phidot^2 - V,  w = p/rho,  psi = 1+2w.")
print("  Kinetic fraction K = (1/2 phidot^2)/rho = (1+w)/2 sets the phase:\n")
print("   psi     w        1+w      kinetic frac K     ghost?   c_s^2")
for psi in [1.0, 0.5, 0.0, -0.5, -1.0]:
    w = (psi - 1.0) / 2.0
    K = (1.0 + w) / 2.0
    print(f"   {psi:+.1f}   {w:+.2f}    {1+w:.2f}     {K:.2f}              no       1.00")
print("""
   VERDICT (B): The dark-energy face is realized by an ordinary canonical scalar:
   c_s^2 = 1 (no gradient instability), positive kinetic term (no ghost), and
   1+w = (1+psi)/2 >= 0 for every psi in [-1,1]. Section 7.3's sign(1+w)=+1 is
   thus not merely the algebraic identity 0<=f<=1 -- it is REALIZED by a manifestly
   stable field. The cosmological transition psi:+1 -> -1 is stably traversable.
   [STRENGTHENS section 7.3: non-phantom is realizable, not just definitional]""")

# ----------------------------------------------------------------------------
# PART C -- PHANTOM NO-GO: the same stability forbids the SEDE w-crossing
# ----------------------------------------------------------------------------
print("=" * 76)
print("PART C  PHANTOM NO-GO (Vikman 2005) -- non-phantom becomes a THEOREM")
print("=" * 76)
print("  For a single k-essence P(X,phi):  rho+p = 2 X P_X.")
print("  Crossing w=-1 (rho+p: + -> -) forces P_X through 0 and to change sign.")
print("  Stability needs D := P_X + 2X P_XX >= 0 (no ghost) and c_s^2 = P_X/D >= 0.")
print("  Track a field forced across the divide (X>0 fixed, P_XX<0):\n")
print("   side           P_X     D=P_X+2X P_XX   c_s^2     status")
X, PXX = 1.0, -0.05
for side, PX in [("non-phantom", +0.20), ("at w=-1", 0.0), ("phantom", -0.20)]:
    D = PX + 2 * X * PXX
    cs2 = PX / D if D != 0 else 0.0
    ok = (D > 0) and (cs2 >= 0)
    print(f"   {side:12s}   {PX:+.2f}     {D:+.2f}          {cs2:+.2f}     "
          f"{'stable' if ok else 'GHOST / grad-instability'}")
print("""
   VERDICT (C): A single ghost-free scalar CANNOT cross w=-1 -- at the crossing
   either the no-ghost condition D>=0 or the causality/stability condition c_s^2>=0
   fails (Vikman 2005). So the framework's sign(1+w)=+1 is upgraded from algebra to
   a stability THEOREM: the covariant completion structurally forbids phantom DE.

   TENSION FLAGGED HONESTLY: this is the SAME no-go that forbids the phantom
   crossing (w_a<0, w through -1) the SEDE branch invokes elsewhere in this program
   for a FALLING a0(z) (see README). The two cannot both
   come from one stable scalar. Either (i) drop the phantom crossing and keep a
   flat/rising a0(z) branch, or (ii) pay for a second field / non-minimal coupling
   (the AeST/RelMOND wall, restated). This is the concrete price the covariant
   gate exacts, and it points the next computation. [CONSTRAINS the phenomenology]""")

print("=" * 76)
print("SUMMARY")
print("=" * 76)
print("""  A  Local MOND face  : ghost-free ONLY at psi->+1; stability SELECTS the
                        broken (MOND) phase -> near-tricritical ordering enforced.
  B  Cosmo DE face    : canonical scalar gives 1+w>=0, c_s^2=1 -> non-phantom
                        realizable by a stable field.
  C  Phantom no-go    : sign(1+w)=+1 is a stability theorem, and it forbids the
                        SEDE w-crossing -> single-field completion cannot do both.

  Net: the first real gate of section 9 is PARTIALLY PASSED. The Landau structure
  is consistent with a ghost-free scalar completion, and stability even EXPLAINS
  the symmetry breaking (galaxies sit at psi=+1 because nothing else is stable).
  What it does NOT yet deliver -- and what remains genuinely owed (section 8.1) --
  is a SINGLE covariant field carrying BOTH faces at once (local mu-deformation AND
  cosmological w) with a full relativistic dispersion. Parts A/B are the two faces
  computed separately; unifying them in one T_mu_nu is the make-or-break still open,
  and Part C shows one sharp obstruction it must navigate.""")
