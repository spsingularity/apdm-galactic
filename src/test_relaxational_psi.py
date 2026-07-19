"""STEP 3: give the order parameter psi the RIGHT dynamics (relaxational, not conservative)
and test whether 'why now' emerges.

The single-scalar failure taught the lesson: psi=1+2w is an emergent ORDER PARAMETER, so
it does NOT obey a conservative wave equation psiddot+3H psidot+V'=0 (that illegitimate
1/2 psidot^2 kinetic term is exactly what broke the gate). An order parameter obeys
DISSIPATIVE relaxational dynamics (Hohenberg-Halperin 'model A'):

      dpsi/dt = -Gamma * dF/dpsi (+ noise),     F(psi;r) = r psi^2 + psi^6,

first order in time, relaxing toward the instantaneous free-energy minimum. The natural
cosmological dissipation is Hubble: Gamma = gamma H. In e-folds (t -> ln a):

      dpsi/dln a = -gamma (2 r psi + 6 psi^5).

The control parameter r(a) is set by the cosmological dark balance,
      r(a) = k (f(a) - 1/2),   f = rho_DE/(rho_DM+rho_DE),
so r<0 (broken / matter phase, psi != 0) at high z, r=0 (critical) at rho_DM=rho_DE
(z~0.3), r>0 (symmetric, psi->0) in the future.

  PART 1  the correct (relaxational) dynamics vs the wrong (conservative) one.
  PART 2  integrate psi(z); it tracks equilibrium at high z.
  PART 3  CRITICAL SLOWING DOWN near r=0: tau_psi/t_H diverges -> psi FREEZES.
          The freeze-out band (tau_psi > t_H) is z ~ (0.09, 0.54) = NOW.
  PART 4  'why now' verdict: we live in the freeze-out band of the DM<->DE phase
          transition; psi is out of equilibrium today. Honest about what is derived
          (the freeze-out structure/width, out-of-equilibrium-now) vs input (the band
          CENTER z~0.3 = rho_DM=rho_DE = the Lambda value).

Self-contained (numpy+scipy). Ref: P. C. Hohenberg, B. I. Halperin, Rev. Mod. Phys. 49,
435 (1977) -- dynamic critical phenomena, models A-J.
"""
import numpy as np
from scipy.integrate import solve_ivp

Om, OL = 0.315, 0.685
k, gamma = 1.0, 1.0


def f_de(z):
    x = (1 + z) ** 3
    return OL / (OL + Om * x)


def r_of_z(z):
    return k * (f_de(z) - 0.5)


def psi_eq(z):
    r = r_of_z(z)
    return (-r / 3) ** 0.25 if r < 0 else 0.0        # broken min if r<0, else symmetric


# ============================================================================
print("=" * 78)
print("PART 1  The correct dynamics: relaxational (model A), NOT conservative")
print("=" * 78)
print("""  WRONG (a fundamental field):   psiddot + 3H psidot + V'(psi) = 0   (2nd order, inertial)
  RIGHT (an order parameter):    dpsi/dt = -Gamma dF/dpsi           (1st order, dissipative)
  An average has no inertia; its effective dynamics is relaxational (Hohenberg-Halperin
  model A). Cosmological dissipation is Hubble: Gamma = gamma H. This is why 'writing
  V(psi) into an action' was illegitimate -- the physical dynamics is a gradient flow.""")

# ============================================================================
print("=" * 78)
print("PART 2/3  Integrate psi(z); tracking at high z, critical slowing down near z~0.3")
print("=" * 78)

def rhs(lna, y):
    z = np.exp(-lna) - 1
    r = r_of_z(z)
    return [-gamma * (2 * r * y[0] + 6 * y[0] ** 5)]

z0 = 3.0
sol = solve_ivp(rhs, [np.log(1 / (1 + z0)), 0.0], [psi_eq(z0)], dense_output=True,
                rtol=1e-8, atol=1e-10, max_step=0.02)

print("   z     r(z)      psi_eq    psi_dyn   tau_psi/t_H   state")
for z in [3, 2, 1, 0.5, 0.3, 0.1, 0.0]:
    pd = float(sol.sol(np.log(1 / (1 + z)))[0])
    r = r_of_z(z)
    Fpp = 8 * abs(r) if r < 0 else 2 * r             # curvature at the equilibrium min
    tau = 1.0 / (gamma * Fpp) if Fpp > 0 else np.inf  # tau_psi/t_H = 1/(gamma F'')
    state = "FROZEN (tau>t_H)" if tau > 1 else "tracks equil."
    print(f"   {z:.1f}   {r:+.3f}    {psi_eq(z):.3f}     {pd:.3f}     {tau:8.2f}     {state}")

rc = 1.0 / (8 * gamma)                                # |r| where tau_psi = t_H
z_of_f = lambda fv: (OL * (1 / fv - 1) / Om) ** (1 / 3) - 1
print(f"\n   critical crossing r=0 (rho_DM=rho_DE):  z = {z_of_f(0.5):.2f}")
print(f"   freeze-out band tau_psi>t_H (|r|<1/(8gamma)={rc:.3f}):  "
      f"z in ({z_of_f(0.5 + rc/k):.2f}, {z_of_f(0.5 - rc/k):.2f})  = NOW (last ~5 Gyr)")

# ============================================================================
print("=" * 78)
print("PART 4  'Why now' verdict")
print("=" * 78)
print("""  The order parameter TRACKS equilibrium at high z (tau_psi << t_H), then FREEZES near
  the rho_DM=rho_DE crossing where critical slowing down sends tau_psi -> infinity: with
  Hubble dissipation, tau_psi ~ t_H exactly around z~0.3, so psi cannot follow equilibrium
  melting to 0 and stays stuck at psi~0.5 TODAY -- the dark sector is OUT OF EQUILIBRIUM
  now, mid-transition.

  So 'why now' becomes: we live in the FREEZE-OUT BAND of the DM<->DE order-parameter
  phase transition (z ~ 0.1-0.5), a dynamically special epoch, rather than a bare
  numerical coincidence. This also splits the two clocks of the earlier downstream Task C:
    - LOCAL psi (galaxies, high density -> r<0 always) condensed at STRUCTURE FORMATION
      (z~1-3, where the background still tracked equilibrium, psi_eq~0.6) and is pinned at
      the broken MOND phase (stability-selected psi=+1, prior turn);
    - COSMOLOGICAL psi (background) is frozen mid-transition NOW.

  HONEST BOUNDARY: the dynamics DERIVES the freeze-out structure -- that the transition is
  critically slowed, out-of-equilibrium, and smeared over Delta z ~ 0.5 around the present
  -- but the band CENTER (z~0.3) is still set by rho_DM=rho_DE, i.e. by the Lambda value,
  which remains input (NOT the retracted 'Lambda^4 = condensation free energy'; only the
  TIMING/tracking is claimed, never the magnitude). So 'why now' is upgraded from a bare
  coincidence to 'we observe during the critically-slowed freeze-out of the dark phase
  transition' -- a partial, honest dynamical explanation, and a soft prediction that the
  dark sector is out of equilibrium (psi frozen ~0.5) today. [STEP 3: partial derivation]""")
