#!/usr/bin/env python3
"""
Thermal-bath consistency of the D-3' KMS-tilt modified-inertia MOND.
If inertia is the reaction to the Deser-Levin/Unruh thermal tilt, the accelerated
system must sample a thermal bath at T_U(a). Count the bath modes inside a galaxy.
(Mirrors APDM's own reframe_test_results falsification of the horizon-thermal route.)
"""
import numpy as np
hbar=1.055e-34; c=3e8; kB=1.38e-23; a0=1.2e-10; H0=2.2e-18
T_U=hbar*a0/(2*np.pi*c*kB); T_dS=hbar*H0/(2*np.pi*kB)
lam_th=hbar*c/(kB*T_U); R_gal=10*3.086e19
Nmodes=(R_gal/lam_th)**3
print(f"T_U(a0)={T_U:.2e} K  T_dS={T_dS:.2e} K  lambda_th={lam_th:.2e} m (~horizon)")
print(f"# thermal bath modes in a 10 kpc galaxy = {Nmodes:.2e}  << 1")
print("=> no thermal bath in a galaxy; the KMS-tilt inertia survives ONLY as a VACUUM")
print("   (not thermal-equilibrium) effect -- the unresolved D-4' covariant-action question.")
