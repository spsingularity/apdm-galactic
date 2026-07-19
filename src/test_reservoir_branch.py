"""KILL TEST 1(i) of the 'reservoir branch' -- a0 as a chemically pinned STATE.

Hypothesis under test (proposed as the salvage of the epsilon/boundary-condition idea in
with the reservoir made physical): the condensate edge
is in CHEMICAL (particle-exchange, not thermal) equilibrium with the uncondensed halo envelope,
whose density is pinned by spherical collapse to the universal rho_vir ~ Delta_vir*rho_crit(z).
Then mu_hat at the edge is universal and the acceleration at which the phonon force leaves the
MOND regime is pinned cosmologically -> a0 would be a state, not a Lagrangian coupling.

The one acceleration the reservoir can pin: in BK's MOND regime (BK 2015; a_phi=sqrt(a0*g_b),
a0=alpha^3 Lambda^2/M_Pl) the phonon gradient satisfies (grad theta)^2 = alpha*M_Pl*g_b.
The MOND form holds while the gradient term dominates X; the crossover is at
    (grad theta)^2 / 2m = m*mu_hat   =>   g_c = 2*m*mu_hat/(alpha*M_Pl).
With the dilute-BEC EOS mu_hat = u0*rho/m the particle mass CANCELS:
    g_c = 2*u0*rho_edge/(alpha*M_Pl).
Only u0 (capped by sigma/m bounds) and rho_edge (the hypothesis) enter.

PASS bar (pre-registered): g_c within ~1 order of magnitude of a0_obs=1.2e-10 m/s^2 at the
repo's derived parameter point, without dialing. Verdict decided by the numbers below.
"""
import numpy as np

# ---- units (consistent with test_geometric_anomaly.py) ----
hbarc = 1.973269804e-5                      # eV*cm
cm = 1.0/hbarc                              # eV^-1
s = 2.99792458e10*cm
Mpc = 3.0856775814913673e21*1e3*cm
Mpl = 2.435e27                              # reduced Planck mass [eV]
H0 = 67.0*(1e5*cm)/s/Mpc                    # [eV]
m_len = 1.0/1.973269804e-7; s_time = 2.99792458e8*m_len; ms2 = m_len/s_time**2
a0_obs = 1.2e-10*ms2                        # observed a0 [eV]
eV4_per_GeVcc = 1e9*hbarc**3                # 1 GeV/cc in eV^4 (1 cm^-3 = hbarc^3 eV^3)
cm2g_per_eV3 = (1.973269804e-5)**2/1.78266192e-33   # sigma/m: 1 eV^-3 -> cm^2/g

# ---- the repo's derived parameter point (paper sec.3-4; compute_solidification.py) ----
m      = 1.0                                # DM boson mass [eV]
alpha_chi = 1.5e-19; m_phi = 0.33e-3        # dark fine structure; mediator [eV]
u0     = 4*np.pi*alpha_chi/m_phi**2         # bare contact repulsion [eV^-2] (paper sec.4)
alpha_BK = 2.5                              # BK's O(1) baryon-phonon coupling (BK 2015 fiducial)

rho_crit = 3*Mpl**2*H0**2                   # [eV^4]
# C=1 onset density: n_c = zeta(3/2)/lambda_dB^3, lambda_dB = 2*pi/(m*v)
v150 = 150*(1e5*cm)/s                       # 150 km/s in units of c
lam_dB = 2*np.pi/(m*v150)                   # [eV^-1]
rho_C1 = m*2.612/lam_dB**3                  # [eV^4]
RESERVOIRS = [
    ("virial envelope 200*rho_crit (the hypothesis)",  200*rho_crit),
    ("C=1 phase boundary at v=150 km/s",               rho_C1),
    ("central galactic 0.3 GeV/cc (generous ceiling)", 0.3*eV4_per_GeVcc),
]

print("="*78)
print("KILL TEST 1(i): reservoir-pinned crossover g_c = 2*u0*rho_edge/(alpha*M_Pl) vs a0")
print("="*78)
print(f"  parameter point: m={m} eV, u0=4*pi*alpha_chi/m_phi^2={u0:.2e} eV^-2, alpha={alpha_BK}")
print(f"  a0_obs = {a0_obs:.3e} eV = 1.2e-10 m/s^2;  rho_crit = {rho_crit:.2e} eV^4")
mu_gal = u0*(0.3*eV4_per_GeVcc)/m
C_gal  = (0.3*eV4_per_GeVcc)/rho_C1
print(f"  [sanity vs paper sec.3-4] mu_int(0.3 GeV/cc)={mu_gal:.1e} eV (paper: ~4e-17);"
      f"  C_gal={C_gal:.1e} (paper: 1.8e6)")
print()
print(f"  {'reservoir':<48}{'rho_edge [eV^4]':>16}{'mu_hat [eV]':>13}{'g_c [m/s^2]':>13}{'g_c/a0':>10}")
for name, rho in RESERVOIRS:
    mu_hat = u0*rho/m
    g_c = 2*m*mu_hat/(alpha_BK*Mpl)         # [eV]
    print(f"  {name:<48}{rho:>16.2e}{mu_hat:>13.2e}{g_c/ms2:>13.2e}{g_c/a0_obs:>10.1e}")

print()
print("="*78)
print("RESCUE ANALYSIS: can ANY parameter point close the gap?  (g_c depends only on u0*rho)")
print("="*78)
rho_res = 200*rho_crit
u0_need = a0_obs*alpha_BK*Mpl/(2*rho_res)
print(f"  u0 needed for g_c=a0 at the virial reservoir: {u0_need:.1e} eV^-2"
      f"  (= {u0_need/u0:.1e} x the derived point)")
sigma_m_now  = u0**2*m/(2*np.pi)*cm2g_per_eV3
sigma_m_need = u0_need**2*m/(2*np.pi)*cm2g_per_eV3
print(f"  implied low-v contact sigma/m: now {sigma_m_now:.1f} cm^2/g -> needed {sigma_m_need:.1e} cm^2/g")
print(f"  (SIDM fits want ~1-5; Bullet needs <~0.5 at cluster v) -> overshoot by"
      f" ~{np.log10(sigma_m_need/5):.0f} orders.")
u0_max = np.sqrt(2*np.pi*(5/cm2g_per_eV3)/m)          # sigma/m(low v) <= 5 cm^2/g
g_c_max = 2*u0_max*rho_res/(alpha_BK*Mpl)
print(f"  sigma/m ceiling (<=5 cm^2/g) caps u0 <= {u0_max:.1e} eV^-2 (~= the derived point!)")
print(f"  -> max attainable g_c = {g_c_max/ms2:.1e} m/s^2 = {g_c_max/a0_obs:.1e} * a0,"
      f"  for ANY m~eV point.")
print(f"  m-scaling: u0_max ~ m^-1/2 -> g_c_max ~ m^-1/2: even m=1e-22 eV (fuzzy DM, a different")
print(f"  theory) gains only ~1e11 -> still ~2 orders short. No corner of parameter space closes it.")

print()
print("="*78)
print("VERDICT")
print("="*78)
gap = a0_obs/(2*u0*rho_res/(alpha_BK*Mpl))
print(f"  g_c(virial reservoir) sits {np.log10(gap):.0f} ORDERS OF MAGNITUDE below a0_obs.")
print("  The chemical potential at any physical reservoir density is far too small to set the")
print("  RAR knee in the BK EFT; the knee is set by the Lagrangian scale a0=alpha^3*Lambda^2/M_Pl")
print("  (a coupling), and the sigma/m ceiling forbids raising u0 to compensate (would need")
print("  ~26 orders more cross-section). The 'reservoir branch' is FALSIFIED as a mechanism")
print("  for a0. [REAL][HIGH] at this parameter point; robust across m via g_c ~ u0*rho only.")
print("  (Caveat: standard BK MOND-regime estimate; Mistele 1909.05710's finite-mu corrections")
print("  change O(1) factors and regime boundaries, not 13 orders.)")
