import numpy as np
hbarc=1.973269804e-5; cm=1.0/hbarc; s=2.99792458e10*cm
GeV=1e9; kpc=3.0856775814913673e21*cm; Mpc=1e3*kpc
Mpl=2.435e27; zeta32=2.6123753486854883; c_kms=2.99792458e5
K=8.617333262e-5
m=1.0
kms=lambda x:x/c_kms
rho_g=0.30*GeV/cm**3; n0=rho_g/m; v_g=kms(150)
alpha=1.53e-19; mphi=kms(100)*m           # derived SIDM point
u0=4*np.pi*alpha/mphi**2                    # induced (=natural) quartic scale

print("### G3 test: is the 'macroscopic healing length' real? ###")
mu_int=u0*n0                                # interaction chem. potential mu=u n
xi_heal=1.0/np.sqrt(2*m*mu_int)            # true healing length
print(f"  u0 = 4pi.alpha/mphi^2      = {u0:.2e} eV^-2")
print(f"  mu_int = u0 n0             = {mu_int:.2e} eV")
print(f"  xi_heal = 1/sqrt(2 m mu)   = {xi_heal/cm:.2e} cm  = {xi_heal/(kpc):.2e} kpc")
print(f"  T_DM (kinetic)             = {0.5*m*v_g**2:.2e} eV   (mu_int/T_DM = {mu_int/(0.5*m*v_g**2):.1e})")
print(f"  ratio xi_heal / (10 kpc)   = {xi_heal/(10*kpc):.1e}   -> healing length is MICROSCOPIC")
print("  => the 10 kpc is the CONDENSED-REGION size (where C>1), NOT the healing length.")
print("     The earlier 5e-39 'tuning' came from forcing xi_heal=10 kpc. Artifact.\n")

print("### honest fragmentation tuning with a physical sextic v ~ u0^2/m_phi ###")
v_sext=u0**2/mphi                           # dimensional three-body estimate [eV^-5]
frag_window=2*v_sext*n0/u0                   # |u_eff|_max/u0 for no fragmentation
print(f"  v ~ u0^2/mphi              = {v_sext:.2e} eV^-5")
print(f"  clean window 2 v n0/u0     = {frag_window:.1e}   (vs artifact 5e-39)")
print("  -> if one INSISTS on first order, tuning is ~1e-13 (sextic-scale dependent), not 1e-39.\n")

print("### the economical resolution: cluster exclusion is KINETIC, so first order is optional ###")
# thermalization already excludes clusters (computed before): cluster Gamma t_dyn = 0.013 < 1
print("  thermalization split alone: group=11 (>1), cluster=0.013 (<1)  [from compute_solidification.py]")
print("  => continuous (XY) branch u0^bare>induced works: no droplets, LSBs clean, clusters kinetic.")
print("     G3 fine-tuning ELIMINATED on the continuous branch.\n")

print("### G1: which horizon quantity sets a_0 ?  (decides the sign of a_0(z)) ###")
H0=67.0*(1e5*cm)/s/Mpc
T_AH=H0/(2*np.pi)                           # Gibbons-Hawking horizon temperature
OmDE=0.689; rho_SEDE=3*Mpl**2*H0**2*OmDE; Lam=rho_SEDE**0.25
# a_0 in eV
m_len=1.0/1.973269804e-7; s_time=2.99792458e8*m_len; ms2=m_len/s_time**2
a0=1.2e-10*ms2
print(f"  a_0(Milgrom)               = {a0:.3e} eV")
print(f"  branch (a) T_AH=H0/2pi     = {T_AH:.3e} eV   (a0/T_AH={a0/T_AH:.2f})")
print(f"  branch (b) sqrt(rhoSEDE)/Mpl = {Lam**2/Mpl:.3e} eV   (a0 = this /{Lam**2/Mpl/a0:.1f} ~ 2pi)")
print(f"  condensate gap mu_int      = {mu_int:.2e} eV ;  mu_int/T_AH = {mu_int/T_AH:.1e}")
print("  => a_0 is NOT the condensate thermal gap (17 orders off T_AH): a_0 is the EFT scale Lam^2/Mpl.")
print("     BK derive a_0~Lam^2/Mpl (phonon-baryon, gravitational) -> DENSITY branch (b), not (a).")
print("     Today (b)~2pi.T_AH by ISW coincidence rho_DE~Mpl^2 H0^2; they DIVERGE in z:")
# z-evolution of both branches using SEDE
import sys
try:
    sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))  # vendored sede/ beside this script
    from sede.friedmann import E_SEDE_lambda
    Om,gam,Or=0.311,1.4964,9.0e-5
    E=lambda z:float(E_SEDE_lambda(np.array([float(z)]),Om,gam,0.5,Or)[0])
    rhoDE=lambda z:E(z)**2-Om*(1+z)**3-Or*(1+z)**4
    for z in [0,1,2,3]:
        a_temp=E(z)/E(0)                     # branch (a): a0 ~ H(z)
        a_dens=np.sqrt(rhoDE(z)/rhoDE(0.0))  # branch (b): a0 ~ sqrt(rho_DE)
        print(f"    z={z}: (a) a0/a0_0 = {a_temp:.3f}   (b) a0/a0_0 = {a_dens:.3f}")
    print("  branch (a) RISES steeply (disfavored); branch (b) FALLS mildly (SEDE prediction). ")
except Exception as e:
    print("   [SEDE unavailable]", e)
