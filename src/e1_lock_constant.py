import numpy as np
# natural units, eV
Mpl=2.435e27
# a0 (Milgrom) in eV
m_len=1.0/1.973269804e-7; s_time=2.99792458e8*m_len; ms2=m_len/s_time**2
a0=1.2e-10*ms2
# H0 = 67 km/s/Mpc in eV
Mpc=3.0856775814913673e24*(1.0/1.973269804e-5)   # Mpc in 1/eV
H0=67.0*(1e5*(1.0/1.973269804e-5))/(2.99792458e10*(1.0/1.973269804e-5))/Mpc
OmDE=0.689
rho_SEDE=3*Mpl**2*H0**2*OmDE
T_AH=H0/(2*np.pi)

print("### The observable lock constant  K = a0 * Mpl / sqrt(rho_SEDE) ###")
K_obs=a0*Mpl/np.sqrt(rho_SEDE)
print(f"  a0            = {a0:.4e} eV")
print(f"  H0            = {H0:.4e} eV   T_AH=H0/2pi = {T_AH:.4e} eV")
print(f"  sqrt(rho_SEDE)/Mpl = {np.sqrt(rho_SEDE)/Mpl:.4e} eV")
print(f"  K_observed    = {K_obs:.4f}")

print("\n### Closed-form derivation (kappa=1 convention: Lambda^4 = rho_SEDE) ###")
# a0 = T_AH = H0/2pi  (horizon-thermodynamic normalization) ; rho_SEDE=3 OmDE Mpl^2 H0^2
K_pred=1.0/(2*np.pi*np.sqrt(3*OmDE))
print(f"  K_pred = 1/(2 pi sqrt(3 Omega_DE)) = {K_pred:.4f}")
print(f"  ratio K_obs/K_pred = {K_obs/K_pred:.3f}   (= a0/T_AH, the residual of Milgrom's coincidence)")
print(f"  a0/T_AH = {a0/T_AH:.3f}")

print("\n### Consequences in the kappa=1 convention ###")
Lam=rho_SEDE**0.25
xi_BK=K_pred                      # predicted BK MOND coefficient a0 = xi_BK * Lambda^2/Mpl
print(f"  Lambda = rho_SEDE^1/4 = {Lam*1e3:.3f} meV")
print(f"  predicted BK coefficient  xi_BK = a0 Mpl/Lambda^2 = {K_pred:.4f}  (NOT free: closed form)")
print(f"  a0_pred = xi_BK * Lambda^2/Mpl = {xi_BK*Lam**2/Mpl:.4e} eV  vs a0_obs={a0:.4e} eV")
# BK write a0 = alpha^3 Lambda^2/Mpl (O(1) prefactor set to 1): implied coupling
alpha_BK=K_obs**(1/3.)
print(f"  if BK form a0=alpha^3 Lambda^2/Mpl (kappa=1): alpha = K_obs^(1/3) = {alpha_BK:.3f}  (~1/2)")

print("\n### Falsifiable closed-form relation (no free parameter) ###")
print(f"  2 pi * K_obs * sqrt(3 Omega_DE) = {2*np.pi*K_obs*np.sqrt(3*OmDE):.3f}  (=1 iff a0=cH0/2pi exactly)")
print("  -> ties a0, Omega_DE, H0 with zero free parameters; holds to 16%.")
