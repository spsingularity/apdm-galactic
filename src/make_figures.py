"""Generate three diagnostic a0(z)/thermalization figures into src/figures_out/.

(These are the APDM-side diagnostic figures; the manuscript's own figures live in
paper/figures/ and are built separately. This script writes to its own output dir so
it never touches the manuscript.)

Figure numbers follow the paper's order of appearance:
  Fig 1 (§3)  the thermalization split Gamma_th t_dyn(v) + sigma/m(v).
  Fig 2 (§5)  a0(z)/a0(0): the fixed-sign three-branch discriminant (density-branch BAND
              spanning SEDE's structural->dynamical w, vs constant, vs rate branch).
  Fig 3 (§7)  'why now' as a critically-slowed freeze-out of the order parameter.

The a0(z) panel imports SEDE (../SEDE) live for w(z) (CPL fallback); the others are self-contained.
Run:  python3 experiments/make_figures.py
"""
import os, sys, json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures_out")  # self-contained; does NOT touch paper/figures
os.makedirs(FIGDIR, exist_ok=True)
sys.path.insert(0, HERE)  # vendored sede/ beside this script
Om, OL = 0.315, 0.685
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3,
                     "figure.dpi": 150, "savefig.bbox": "tight"})

def save_both(fig, name):
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIGDIR, f"{name}.{ext}"))
    plt.close(fig)

# ============================================================================
# Fig 2 (paper §5) -- a0(z)/a0(0): the fixed-sign discriminant
# ============================================================================
zg = np.linspace(0, 3, 400)
lna = np.log(1.0 / (1 + zg))

def a0_ratio(wz):                                  # a0 ∝ sqrt(rho_DE),  dln rho/dln a = -3(1+w)
    integ = np.cumsum((-3 * (1 + wz))[::-1] * np.gradient(lna[::-1]))[::-1]
    return np.sqrt(np.exp(integ - integ[0]))

try:
    from sede import theory
    pj = json.load(open(os.path.join(HERE, "predictions.json")))
    w_mild = pj["w0"] + pj["wa"] * zg / (1 + zg)                       # structural CPL
    w_dyn = np.array([float(np.atleast_1d(theory.w_DE_dynamical(np.array([z]), 0.311))[0]) for z in zg])
    tag = "SEDE"
except Exception as e:
    print(f"[SEDE unavailable ({e}); CPL fallback]")
    w_mild = -0.98 - 0.11 * zg / (1 + zg); w_dyn = -1.0 - 0.30 * zg / (1 + zg); tag = "CPL"

mild, steep = a0_ratio(w_mild), a0_ratio(w_dyn)
rate = np.sqrt(Om * (1 + zg) ** 3 + OL)

fig, ax = plt.subplots(figsize=(7, 5))
# current systematics floor: +-0.2 dex around unity (can't distinguish from constant inside)
ax.axhspan(10 ** -0.2, 10 ** 0.2, color="0.85", zorder=0, label="±0.2 dex systematics floor")
ax.fill_between(zg, steep, mild, color="C0", alpha=0.35, zorder=2,
                label="density branch $a_0\\!\\propto\\!\\sqrt{\\rho_{\\rm DE}}$ (this work): FALLING")
ax.plot(zg, mild, color="C0", lw=1.2, ls="--", alpha=0.8)
ax.plot(zg, steep, color="C0", lw=1.2, ls="--", alpha=0.8)
ax.plot(zg, np.ones_like(zg), color="k", lw=2, label="constant $a_0$ (standard MOND / ΛCDM)")
ax.plot(zg, rate, color="C3", lw=2, label="rate branch $a_0\\!\\propto\\! cH(z)$ (Verlinde): RISING")
ax.axvspan(1, 3, color="C2", alpha=0.08, zorder=0)
ax.text(2.0, 0.30, "JWST / Euclid\nhigh-$z$ reach", color="C2", ha="center", fontsize=9)
ax.set_yscale("log")
ax.set_xlabel("redshift $z$"); ax.set_ylabel("$a_0(z)/a_0(0)$")
ax.set_xlim(0, 3); ax.set_ylim(0.28, 5.2)
ax.set_yticks([0.3, 0.5, 0.7, 1, 2, 3, 5]); ax.set_yticklabels(["0.3", "0.5", "0.7", "1", "2", "3", "5"])
ax.set_title(f"The MOND scale vs redshift: the sign is the test  ({tag})")
ax.legend(fontsize=8.5, loc="upper left", framealpha=0.9)
save_both(fig, "fig2_a0_of_z")
print("wrote fig2_a0_of_z .png/.pdf")

# ============================================================================
# Fig 1 (paper §3) -- the thermalization split (the §3 table, made visual)
# ============================================================================
sys_name = ["dwarf", "LSB", "HSB", "group", "cluster"]
v = np.array([30, 70, 150, 500, 1200.0])
Gam = np.array([4.3e6, 8.5e4, 1.9e5, 11.0, 0.013])
sig = np.array([4.6, 3.0, 1.0, 0.038, 0.002])
col = ["C0", "C0", "C0", "C1", "C3"]

fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 4.4))
a1.axhline(1.0, color="k", ls="--", lw=1.5)
a1.text(33, 1.6, "$\\Gamma_{\\rm th}t_{\\rm dyn}=1$  (condensation threshold)", fontsize=8.5)
for i in range(5):
    a1.scatter(v[i], Gam[i], s=90, color=col[i], zorder=3, edgecolor="k", lw=0.5)
    a1.annotate(f"{sys_name[i]}\n({Gam[i]:.2g})", (v[i], Gam[i]),
                textcoords="offset points", xytext=(0, 10 if i != 3 else -28), ha="center", fontsize=8.5)
a1.scatter([], [], color="C0", label="superfluid"); a1.scatter([], [], color="C1", label="crossover (APDM)")
a1.scatter([], [], color="C3", label="normal")
a1.set_xscale("log"); a1.set_yscale("log")
a1.set_xlabel("velocity dispersion $v$ [km/s]"); a1.set_ylabel("$\\Gamma_{\\rm th}\\, t_{\\rm dyn}$")
a1.set_title("Thermalization split: the group-scale crossover"); a1.legend(fontsize=9, loc="upper right")
a1.set_ylim(3e-3, 3e7)

a2.axhline(0.5, color="C3", ls="--", lw=1.5); a2.text(90, 0.6, "Bullet Cluster limit  $\\sigma/m<0.5$", color="C3", fontsize=8.5)
a2.plot(v, sig, "o-", color="0.3", zorder=3)
for i in range(5):
    a2.annotate(sys_name[i], (v[i], sig[i]), textcoords="offset points", xytext=(6, 6), fontsize=8.5)
a2.set_xscale("log"); a2.set_yscale("log")
a2.set_xlabel("velocity $v$ [km/s]"); a2.set_ylabel("$\\sigma/m$ [cm$^2$/g]")
a2.set_title("Velocity-dependent self-interaction (Bullet-safe)")
save_both(fig, "fig1_thermalization_split")
print("wrote fig1_thermalization_split .png/.pdf")

# ============================================================================
# Fig 3 (paper §7) -- 'why now' as a critically-slowed freeze-out
# ============================================================================
k, gamma = 1.0, 1.0
def f_de(z): return OL / (OL + Om * (1 + z) ** 3)
def r_of_z(z): return k * (f_de(z) - 0.5)
def psi_eq(z):
    r = r_of_z(z); return (-r / 3) ** 0.25 if r < 0 else 0.0
def rhs(N, y):
    z = np.exp(-N) - 1; return [-gamma * (2 * r_of_z(z) * y[0] + 6 * y[0] ** 5)]
sol = solve_ivp(rhs, [np.log(1 / 4.0), 0.0], [psi_eq(3.0)], dense_output=True, rtol=1e-8, atol=1e-10, max_step=0.02)
zz = np.linspace(3, 0, 400)
psid = np.array([float(sol.sol(np.log(1 / (1 + z)))[0]) for z in zz])
psie = np.array([psi_eq(z) for z in zz])
tau = np.array([1.0 / (gamma * (8 * abs(r_of_z(z)) if r_of_z(z) < 0 else 2 * r_of_z(z) + 1e-9)) for z in zz])
z_eq = (OL / Om) ** (1 / 3) - 1
zlo = (OL * (1 / (0.5 + 1 / 8) - 1) / Om) ** (1 / 3) - 1
zhi = (OL * (1 / (0.5 - 1 / 8) - 1) / Om) ** (1 / 3) - 1

fig, (b1, b2) = plt.subplots(1, 2, figsize=(11, 4.4))
for ax in (b1, b2):
    ax.axvspan(zlo, zhi, color="C1", alpha=0.18, zorder=0)
    ax.axvline(z_eq, color="C1", ls=":", lw=1.5)
b1.plot(zz, psie, color="C7", lw=2, ls="--", label="equilibrium $\\psi_{\\rm eq}(z)$")
b1.plot(zz, psid, color="C0", lw=2.5, label="dynamical $\\psi(z)$ (model A)")
b1.text(z_eq, 0.66, "  $\\rho_{\\rm DM}=\\rho_{\\rm DE}$\n  ($z\\approx0.3$)", color="C1", fontsize=8.5, va="top")
b1.text((zlo + zhi) / 2, 0.05, "freeze-out\nband = NOW", color="C1", ha="center", fontsize=8.5)
b1.set_xlabel("redshift $z$"); b1.set_ylabel("order parameter $\\psi$")
b1.set_title("The dark sector freezes mid-transition, today"); b1.legend(fontsize=9, loc="lower left")
b1.set_xlim(3, 0)
b2.axhline(1.0, color="k", ls="--", lw=1.5); b2.text(2.9, 1.3, "$\\tau_\\psi=t_H$", fontsize=9)
b2.plot(zz, tau, color="C0", lw=2.5)
b2.text(2.4, 15, "critical slowing down\n($\\tau_\\psi\\!\\to\\!\\infty$ at the crossing)", fontsize=8.5)
b2.set_yscale("log"); b2.set_xlabel("redshift $z$"); b2.set_ylabel("$\\tau_\\psi / t_H$")
b2.set_title("Relaxation time vs Hubble time"); b2.set_xlim(3, 0); b2.set_ylim(0.2, 1e2)
save_both(fig, "fig3_freezeout")
print("wrote fig3_freezeout .png/.pdf")
print("all figures ->", os.path.normpath(FIGDIR))
