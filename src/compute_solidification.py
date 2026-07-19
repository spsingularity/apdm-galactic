import numpy as np
# ---- natural units hbar=c=1, energies in eV ----
hbarc=1.973269804e-5            # eV*cm
cm=1.0/hbarc                    # 1 cm  [1/eV]
s=2.99792458e10*cm             # 1 s   [1/eV]
g_gram=5.6095886e32            # 1 g   [eV]
GeV=1e9; K=8.617333262e-5      # 1 K [eV]
kpc=3.0856775814913673e21*cm; Mpc=1e3*kpc
Mpl=2.435e27                   # reduced Planck mass [eV]
G_N=1.0/(8*np.pi*Mpl**2)       # [1/eV^2]
zeta32=2.6123753486854883
eVm3_to_cm2g=hbarc**2*g_gram   # 1 eV^-3 -> cm^2/g  (=2.185e23)
c_kms=2.99792458e5
def rho_gevcc(x): return x*GeV/cm**3          # GeV/cm^3 -> eV^4
def kms(x): return x/c_kms                    # km/s -> units of c

m=1.0                                          # boson mass [eV]
# ---------- 1. self-consistency of the m=1 eV point ----------
v_g=kms(150); rho_g=rho_gevcc(0.30); n_g=rho_g/m
lam=lambda v:2*np.pi/(m*v)
C=lambda rho,v:(rho/m)*lam(v)**3/zeta32
T_DM=0.5*m*v_g**2
Tc=2*np.pi*(n_g/zeta32)**(2/3)/m
print("### 1. m=1 eV self-consistency")
print(f"  n_gal={n_g:.3e} eV^3   lambda_dB={lam(v_g)/cm:.3f} cm   C_gal={C(rho_g,v_g):.2e}")
print(f"  T_DM={T_DM/K*1e3:.2f} mK   T_c={Tc/K:.2f} K   T_DM/T_c={T_DM/Tc:.1e}")

# ---------- 2. SIDM Born-Yukawa, anchored at sigma/m(150 km/s)=1 cm^2/g ----------
def L(xi): return np.log(1+xi**2)-xi**2/(1+xi**2)
def sigT(v,alpha,mphi): return 8*np.pi*alpha**2/(m**2*v**4)*L(m*v/mphi)   # [eV^-2]
def som_cm2g(v,alpha,mphi): return sigT(v,alpha,mphi)/m*eVm3_to_cm2g
vstar_kms=100.0; mphi=kms(vstar_kms)*m
# solve alpha from sigma/m(150)=1 cm^2/g
tgt=1.0/eVm3_to_cm2g*m                          # sigma target [eV^-2]
alpha=np.sqrt(tgt/(8*np.pi/(m**2*kms(150)**4)*L(m*kms(150)/mphi)))
print("\n### 2. SIDM point (derived: anchor sigma/m(150)=1 cm^2/g)")
print(f"  m_phi={mphi*1e3:.3f} meV (v_*={vstar_kms:.0f} km/s)   alpha_chi={alpha:.3e}   alpha*m/m_phi={alpha*m/mphi:.1e} (Born OK)")
for vv in [30,70,150,500,1200]:
    print(f"    v={vv:5d} km/s   sigma/m={som_cm2g(kms(vv),alpha,mphi):8.3f} cm^2/g")
a_sl=-np.log(som_cm2g(kms(1200),alpha,mphi)/som_cm2g(kms(150),alpha,mphi))/np.log(1200/150)
print(f"  effective slope a(150->1200)={a_sl:.2f}")

# ---------- 3. thermalization split from derived sigma(v) ----------
print("\n### 3. thermalization  Gamma_th*t_dyn = N_occ (sigma/m) rho v * m^? ...")
print("  system   rho[GeV/cc] v[km/s]  Gamma_th*t_dyn  verdict")
for name,rg,vv in [("HSB",0.30,150),("LSB",0.03,70),("dwarf",0.10,30),("group",0.02,500),("cluster",0.005,1200)]:
    rho=rho_gevcc(rg); v=kms(vv); n=rho/m
    Nocc=(2*np.pi)**3*rho/(m**4*v**3)
    som=sigT(v,alpha,mphi)/m                     # sigma/m [eV^-3]
    Gamma=Nocc*(som*m)*n*v                        # = N (sigma) n v ; sigma=(sigma/m)*m
    tdyn=1/np.sqrt(G_N*rho)
    p=Gamma*tdyn
    print(f"  {name:7s}  {rg:8.3f}   {vv:5d}   {p:12.2e}   {'superfluid' if p>1 else 'NORMAL'}")

# ---------- 4. tricritical block; sextic from macroscopic healing length ----------
print("\n### 4. tricritical numbers")
xi_coh=10*kpc; n0=n_g
# effective quartic giving a 10 kpc healing length:  xi=1/sqrt(2 m u n)  => u_eff_edge
u_edge=1.0/(2*m*n0*xi_coh**2)
v_sext=1.0/(4*n0**2*xi_coh**2*m)                 # matched-edge sextic (|u_eff|=2 v n0 at edge)
Kk=1.0/(2*m)
Gi=T_DM*np.sqrt(v_sext)/Kk**1.5
Ncoh=n0*xi_coh**3
Fc_T=0.74/(Gi*0.1**2)
print(f"  u_eff@10kpc edge = {u_edge:.2e} eV^-2   sextic v={v_sext:.2e} eV^-5")
print(f"  Gi_tri={Gi:.2e}   N_in_xi^3={Ncoh:.2e}   F_c/T={Fc_T:.2e}")
# natural Yukawa quartic and the clean-window width
u0=4*np.pi*alpha/mphi**2
print(f"  natural Yukawa quartic u0=4pi*alpha/m_phi^2 = {u0:.2e} eV^-2")
print(f"  clean-window frac  |u_eff|_max/u0 = 2 v n0/u0 = {2*v_sext*n0/u0:.2e}")
print(f"  => equivalently u_eff must sit within {u_edge/u0:.1e} of the tricritical point")

# ---------- 5. horizon numerology ----------
print("\n### 5. horizon / a_0 numerology")
H0=67.0*(1e5*cm)/s/Mpc        # 67 km/s/Mpc  [eV]
OmDE=0.689; rho_SEDE=3*Mpl**2*H0**2*OmDE; Lam=rho_SEDE**0.25
a0=Lam**2/Mpl
m_len=1.0/1.973269804e-7      # 1 m [1/eV]
s_time=2.99792458e8*m_len     # 1 s [1/eV]
ms2=m_len/s_time**2           # 1 m/s^2 [eV]
a0_milg=1.2e-10*ms2
print(f"  Lambda=rho_SEDE^1/4={Lam*1e3:.3f} meV (~m_phi={mphi*1e3:.2f} meV)")
print(f"  a0=Lambda^2/Mpl={a0:.2e} eV   a0(Milgrom)={a0_milg:.2e} eV   ratio={a0/a0_milg:.2f} (~2pi)")

# ---------- 6. corrected a_0(z) from SEDE's actual rho_DE(a) (optional) ----------
try:
    import os, sys
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))  # vendored sede/ beside this script
    from sede.friedmann import E_SEDE_lambda
    Om,gam,Or=0.311,1.4964,9.0e-5
    E=lambda z:E_SEDE_lambda(z,Om,gam,0.5,Or)
    rhoDE=lambda z:E(z)**2-Om*(1+z)**3-Or*(1+z)**4
    zz=np.array([0.,0.5,1.,2.,3.])
    a0z=np.sqrt(rhoDE(zz)/rhoDE(0.))
    print("\n### 6. corrected a_0(z) from SEDE (a0 ~ sqrt(rho_SEDE))")
    for z,a in zip(zz,a0z): print(f"    z={z:.1f}  a0(z)/a0(0)={a:.4f}")
    print("    -> a_0 DECREASES with z (SEDE w crosses -1, w_a<0); opposite sign to a thawing CPL")
except Exception as e:
    print(f"\n### 6. (SEDE not on PYTHONPATH: run with PYTHONPATH=../SEDE)  [{e}]")
