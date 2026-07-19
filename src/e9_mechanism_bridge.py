import numpy as np, sys
sys.path.insert(0, __import__('os').path.dirname(__import__('os').path.abspath(__file__)))  # vendored sede/ beside this script
from sede.friedmann import compute_growth_factor, E_SEDE_lambda
from scipy.integrate import quad
from scipy.special import erfc

Om,h,ns,s8=0.311,0.674,0.965,0.811
rho_m=2.775e11*Om*h**2                     # Msun/Mpc^3 (comoving), rho_crit=2.775e11 h^2
dc=1.686
# --- BBKS transfer + sigma(M) ---
Gam=Om*h
def T(k):                                   # k in 1/Mpc
    q=k/Gam
    return np.log(1+2.34*q)/(2.34*q)*(1+3.89*q+(16.1*q)**2+(5.46*q)**3+(6.71*q)**4)**-0.25
def P(k,A): return A*k**ns*T(k)**2
def sig2(R,A):
    def integ(lk):
        k=np.exp(lk); x=k*R; W=3*(np.sin(x)-x*np.cos(x))/x**3
        return k**3*P(k,A)/(2*np.pi**2)*W**2
    return quad(integ,np.log(1e-4),np.log(1e3),limit=200)[0]
A=1.0; A=s8**2/sig2(8.0/h,1.0)              # normalize sigma8 at 8 Mpc/h
def R_of_M(M): return (3*M/(4*np.pi*rho_m))**(1/3.)
def sigM(M): return np.sqrt(sig2(R_of_M(M),A))
# tabulate sigma(M) z=0
lgM=np.linspace(6,16,60); sig0=np.array([sigM(10**l) for l in lgM])
def sigma0(M): return np.interp(np.log10(M),lgM,sig0)

# --- Sheth-Tormen multiplicity ---
def fST(nu):
    A_,a,p=0.3222,0.707,0.3
    return A_*np.sqrt(2*a/np.pi)*(1+(1/(a*nu**2))**p)*nu*np.exp(-a*nu**2/2)
def dndlnM(M,Dz):
    s=sigma0(M)*Dz
    dl=1e-3; dlnsdlnM=(np.log(sigma0(M*np.exp(dl)))-np.log(sigma0(M*np.exp(-dl))))/(2*dl)
    nu=dc/s
    return rho_m/M*fST(nu)*abs(dlnsdlnM)     # dn/dlnM  (Msun^-1 Mpc^-3 * ... )

Ms=np.logspace(7,15.5,80)
def mass_weighted(Dz,power=1.0):
    integ=np.array([dndlnM(M,Dz)*M**power for M in Ms])
    return np.trapezoid(integ,np.log(Ms))
def Fcoll_PS(Dz,Mmin=1e8):                   # Press-Schechter collapsed fraction
    return erfc(dc/(np.sqrt(2)*sigma0(Mmin)*Dz))

# --- growth + curves ---
zs=np.array([0,0.3,0.5,0.7,1.0,1.5,2.0,3.0])
D=compute_growth_factor(zs,Om); D=D/D[0]     # normalized growth, D(0)=1
gam=1.50
fsat=(1-np.exp(-gam*D**2))/(1-np.exp(-gam))  # SEDE f_sat
# APDM condensate: mass-weighted (S~1) and entropy-weighted (S(M)~M^p, p=1.2)
fcond_mass=np.array([mass_weighted(d,1.0) for d in D]); fcond_mass/=fcond_mass[0]
fcond_ent =np.array([mass_weighted(d,1.2) for d in D]); fcond_ent/=fcond_ent[0]
Fps=np.array([Fcoll_PS(d) for d in D]); Fps/=Fps[0]

print(" z    D(z)   f_sat(SEDE)  f_cond(mass)  f_cond(entropy p=1.2)  F_coll(PS)")
for i,z in enumerate(zs):
    print(f" {z:4.1f} {D[i]:.3f}   {fsat[i]:.3f}       {fcond_mass[i]:.3f}         {fcond_ent[i]:.3f}              {Fps[i]:.3f}")
# quality of match
def chi(a,b): return np.sqrt(np.mean((a-b)**2))
print(f"\n RMS(f_cond_mass - f_sat)    = {chi(fcond_mass,fsat):.3f}")
print(f" RMS(f_cond_entropy - f_sat) = {chi(fcond_ent,fsat):.3f}")
print(f" RMS(F_coll_PS - f_sat)      = {chi(Fps,fsat):.3f}")
# where do they diverge?
print(f"\n ratio f_cond_ent/f_sat at z=1,2,3: {fcond_ent[4]/fsat[4]:.2f}, {fcond_ent[6]/fsat[6]:.2f}, {fcond_ent[7]/fsat[7]:.2f}")
