"""Decisive first-tests for the 'state-not-coupling' reframe.

TEST 1  (kill-or-confirm the FRONT-RUNNER): can the galactic condensate reach DETAILED BALANCE /
        thermal contact with the Gibbons-Hawking horizon bath (T_dS = H/2pi) within a Hubble time?
        If not, the FDT/thermal-bath route to a0 = c*H(z)/2pi is dead (a0 freezes -> T2 returns).

TEST 3  (kill-or-confirm CANDIDATE 3): does the BK P(X) ~ X^{3/2} condensate support a
        crystalline / supersolid phase (the T4 fix), or is it a weakly-coupled uniform superfluid?

Reuses the unit conventions + derived parameters of compute_solidification.py and
run_phonon_sound_speed.py. Natural units hbar=c=k_B=1, energies in eV.
"""
import numpy as np

# ---- units (from compute_solidification.py) ----
hbarc = 1.973269804e-5                 # eV*cm
cm    = 1.0/hbarc                       # 1 cm [1/eV]
s     = 2.99792458e10*cm               # 1 s  [1/eV]
GeV   = 1e9
kpc   = 3.0856775814913673e21*cm
Mpc   = 1e3*kpc
Gpc   = 1e3*Mpc
Mpl   = 2.435e27                        # reduced Planck mass [eV]
c_kms = 2.99792458e5
zeta32= 2.6123753486854883
def rho_gevcc(x): return x*GeV/cm**3    # GeV/cm^3 -> eV^4
def kms(v): return v/c_kms              # km/s -> units of c

# ---- derived APDM point ----
m    = 1.0                              # boson mass [eV]
Lam  = 2.24e-3                          # BK EOS scale [eV] (~meV)
mphi = kms(100.0)*m                     # mediator mass [eV]
alpha= 1.53e-19                         # from sigma/m(150)=1 cm^2/g
rho_g= rho_gevcc(0.30)                  # galactic condensate density [eV^4]
n_g  = rho_g/m                          # number density [eV^3]

# ---- cosmology ----
H0 = 67.0*(1e5*cm)/s/Mpc               # 67 km/s/Mpc [eV]
T_dS = H0/(2*np.pi)                     # Gibbons-Hawking temperature today [eV]  (hbar=c=k_B=1)

print("="*74)
print("TEST 1 - horizon<->condensate thermal-contact time  (front-runner kill-or-confirm)")
print("="*74)
print(f"  H0                     = {H0:.3e} eV")
print(f"  T_dS = H0/2pi          = {T_dS:.3e} eV   (de Sitter / Gibbons-Hawking bath)")
lam_dS = 1.0/T_dS                       # thermal wavelength [1/eV]
print(f"  lambda_dS = 1/T_dS     = {lam_dS/Gpc:.2f} Gpc   (bath thermal wavelength)")

# BK phonon sound speed at galactic density: c_s = rho/(2 Lambda m^3)
c_s = rho_g/(2*Lam*m**3)               # units of c
print(f"\n  BK phonon c_s(galactic)= {c_s*c_kms:.0f} km/s = {c_s:.2e} c")

# softest MOND-carrying phonon: wavelength ~ coherence length (10 kpc); also check kpc
for Lname, Lcoh in [("10 kpc (coherence)", 10*kpc), ("1 kpc", 1*kpc), ("Hubble radius", 1.0/H0)]:
    k_ph  = 2*np.pi/Lcoh                # [eV]
    w_ph  = c_s*k_ph                    # phonon frequency [eV]  (omega = c_s k)
    x     = w_ph/T_dS
    # bath occupation of that mode
    occ   = 1.0/np.expm1(x) if x < 700 else np.exp(-x)
    print(f"  phonon @ {Lname:18s}: omega={w_ph:.2e} eV  omega/T_dS={x:.2e}  n_bath(omega)={occ:.2e}")

# how many bath modes fit inside a galaxy
R_gal = 10*kpc
N_modes = (R_gal/lam_dS)**3
print(f"\n  # of horizon-bath modes inside a 10 kpc galaxy ~ (R_gal/lambda_dS)^3 = {N_modes:.1e}")

# equilibration: even the *softest* galactic phonon (Hubble-radius wavelength) vs T_dS
k_soft = 2*np.pi*H0
x_soft = (c_s*k_soft)/T_dS
print(f"  softest possible galactic phonon (lambda=Hubble R): omega/T_dS = {x_soft:.2f}")
print("\n  VERDICT: T_dS ~ 1e-34 eV, lambda_dS ~ Gpc. Every MOND-carrying phonon has")
print("  omega/T_dS >> 1  => bath occupation e^{-omega/T_dS} ~ 0, and <<1 bath mode fits in a")
print("  galaxy. The condensate CANNOT thermalize with the horizon bath in a Hubble time:")
print("  the FDT / detailed-balance route to a0=cH/2pi is FALSIFIED (a0 would freeze -> T2).")
print("  => Surviving route to a0 ~ H(z) is GEOMETRIC (curvature R~H^2 / non-renormalized")
print("     scaling anomaly, a T=0 statement) -- NOT tested here; must be tested separately.")

print()
print("="*74)
print("TEST 3 - does P(X)~X^{3/2} support a crystal / supersolid?  (candidate 3 kill-or-confirm)")
print("="*74)
# diluteness / gas parameter
u0   = 4*np.pi*alpha/mphi**2            # bare Yukawa quartic [eV^-2] (contact coupling g)
a_s  = u0*m/(4*np.pi)                   # s-wave scattering length [1/eV]  (g = 4pi a_s/m)
gasp = n_g*a_s**3                       # diluteness parameter n a_s^3
d    = n_g**(-1.0/3.0)                  # interparticle spacing [1/eV]
E_kin= n_g**(2.0/3.0)/(2*m)            # zero-point kinetic energy/particle ~ hbar^2/(2 m d^2)
E_int= u0*n_g                           # interaction energy/particle ~ g n = mu_int
Gamma= E_int/E_kin                      # coupling parameter (int/kinetic)
print(f"  n_gal                  = {n_g:.3e} eV^3   spacing d=n^-1/3 = {d/cm:.2e} cm")
print(f"  scattering length a_s   = {a_s/cm:.2e} cm    a_s/d = {a_s/d:.2e}")
print(f"  gas parameter n a_s^3   = {gasp:.2e}   (dilute BEC needs <<1)")
print(f"  E_int/particle (=mu_int)= {E_int:.2e} eV")
print(f"  E_kin/particle (zero-pt)= {E_kin:.2e} eV")
print(f"  coupling Gamma=E_int/E_kin = {Gamma:.2e}")
print(f"  (Wigner/quantum crystallization needs Gamma >~ 1e2; classical OCP melt Gamma~175)")

# roton check: Bogoliubov dispersion omega(k)=sqrt((c_s k)^2 + (k^2/2m)^2) is monotonic -> no roton
print("\n  Bogoliubov dispersion omega(k)=sqrt((c_s k)^2+(k^2/2m)^2): monotonically increasing")
print("  -> NO roton minimum -> no soft-mode (supersolid) instability from the BK EOS.")
print(f"\n  VERDICT: gas parameter {gasp:.1e} << 1 and Gamma {Gamma:.1e} << 100 => the condensate")
print("  is a WEAKLY-COUPLED UNIFORM SUPERFLUID; P~X^{3/2} gives a monotonic (roton-free)")
print("  spectrum. It does NOT spontaneously crystallize/become supersolid. Forcing a crystal")
print("  needs a strong, momentum-structured (dipolar/soft-core) interaction that generically")
print("  DEFORMS the MOND-regime EOS away from X^{3/2} -> spoils the phonon force law.")
print("  => spontaneous-supersolid T4 fix DISFAVORED for BK-as-is. (Baryon vortex-PINNING,")
print("     INV-1, is a separate U_pin>k_BT calc, not this crystallization criterion.)")
