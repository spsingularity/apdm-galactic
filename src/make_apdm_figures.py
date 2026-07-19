"""Generate the two figures for the rewritten APDM paper, DESI-anchored and self-contained.

  Fig 1 (§3)  a0(z)/a0(0): the fixed-sign three-branch discriminant. The density branch is a BAND
              built directly from the DESI DR2 w0waCDM posteriors (Pantheon+, DES-Y5, Union3), with
              the SEDE model drawn as one example curve. Rate branch a0∝cH and constant a0 shown for
              contrast, against the ±0.2 dex systematics floor.
  Fig 2 (§4)  the preregistered forecast: dlog10(sigma) at fixed M_b vs observation redshift z_obs
              for the three branches (density anchored to the DESI-preferred w), and the
              separability-vs-N curve.

Everything is anchored to the DESI DR2 equation of state, so the discriminant is data-driven and
model-independent; SEDE (a sibling model) appears only as one example of an evolving-w background.
Run:  python3 experiments/make_apdm_figures.py            # regenerate the figure PDFs
      python3 experiments/make_apdm_figures.py --no-save  # print the numbers only (no file writes)
"""
import os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
FIGDIR = os.path.join(HERE, "figures_out")  # self-contained; does NOT touch paper/figures
# --no-save (or APDM_NO_SAVE=1) prints every quantitative result but does NOT
# write the figure files, so reproduce.py can verify the numbers without
# writing anything.
NO_SAVE = ("--no-save" in sys.argv) or (os.environ.get("APDM_NO_SAVE") == "1")
if not NO_SAVE:
    os.makedirs(FIGDIR, exist_ok=True)

def save_fig(fig, stem):
    if NO_SAVE:
        return
    for ext in ("png", "pdf"):
        fig.savefig(os.path.join(FIGDIR, f"{stem}.{ext}"))
plt.rcParams.update({"font.size": 11, "axes.grid": True, "grid.alpha": 0.3,
                     "figure.dpi": 150, "savefig.bbox": "tight"})
Om, OL = 0.31, 0.69

# ---- equations of state --------------------------------------------------
# DESI DR2 w0waCDM central values (DESI+CMB + each SN sample), arXiv:2503.14738.
DESI = {"Pantheon+": (-0.838, -0.62), "DES-Y5": (-0.752, -0.86), "Union3": (-0.667, -1.09)}
DESI_REP = DESI["DES-Y5"]                      # representative "DESI-preferred" central (4.2σ)
SEDE = (-0.984, -0.109)                         # sibling model, one example curve

def rho_ratio(z, w0, wa):                       # rho_DE(z)/rho_DE(0) for CPL w(a)=w0+wa(1-a)
    return (1 + z) ** (3 * (1 + w0 + wa)) * np.exp(-3 * wa * z / (1 + z))

def a0_density(z, w0, wa): return np.sqrt(rho_ratio(z, w0, wa))     # density branch a0∝√rho_DE
def a0_rate(z):            return np.sqrt(Om * (1 + z) ** 3 + OL)   # rate branch a0∝cH=E(z)

# ==========================================================================
# Figure 1 -- the fixed-sign discriminant, density band from the DESI posteriors
# ==========================================================================
zg = np.linspace(0, 3, 400)
dens = np.array([a0_density(zg, *DESI[k]) for k in DESI])
lo, hi = dens.min(0), dens.max(0)

fig, ax = plt.subplots(figsize=(7, 5))
ax.axhspan(10 ** -0.2, 10 ** 0.2, color="0.85", zorder=0, label="±0.2 dex systematics floor")
ax.fill_between(zg, lo, hi, color="C0", alpha=0.35, zorder=2,
                label=r"density $a_0\!\propto\!\sqrt{\rho_{\rm DE}}$ (DESI DR2 band)")
ax.plot(zg, a0_density(zg, *SEDE), color="C0", lw=1.3, ls=":", zorder=3, label="  — SEDE example")
ax.plot(zg, np.ones_like(zg), color="k", lw=2, label=r"constant $a_0$ (standard MOND / ΛCDM)")
ax.plot(zg, a0_rate(zg), color="C3", lw=2, label=r"rate $a_0\!\propto\! cH(z)$ (Verlinde): RISES")
ax.axvspan(1, 3, color="C2", alpha=0.08, zorder=0)
ax.text(2.0, 0.62, "JWST / Euclid\nhigh-$z$ reach", color="C2", ha="center", fontsize=9)
ax.set_yscale("log"); ax.set_xlim(0, 3); ax.set_ylim(0.6, 5.2)
ax.set_yticks([0.7, 1, 2, 3, 5]); ax.set_yticklabels(["0.7", "1", "2", "3", "5"])
ax.set_xlabel("redshift $z$"); ax.set_ylabel("$a_0(z)/a_0(0)$")
ax.set_title("The MOND scale vs redshift: the sign is the test")
ax.legend(fontsize=8.5, loc="upper left", framealpha=0.9)
save_fig(fig, "fig_a0_of_z")
plt.close(fig)

# ==========================================================================
# Figure 2 -- the preregistered forecast vs z_obs, density anchored to the DESI-preferred w
# ==========================================================================
y = lambda z, w0, wa: 0.25 * np.log10(a0_density(z, w0, wa))       # dlog10(sigma), density branch
yR = lambda z: 0.25 * np.log10(a0_rate(z))                        # dlog10(sigma), rate branch
sigma_dex = 0.05
rng = np.random.default_rng(0)
zgrid = np.linspace(0.0, 3.0, 61)
yd = np.interp                                                     # (linear on grid)
gD = np.array([y(z, *DESI_REP) for z in zgrid])
gR = np.array([yR(z) for z in zgrid])
gC = np.zeros_like(zgrid)

def forecast(N, grid, ntrial=4000):
    zt = rng.uniform(0.3, 2.5, size=(ntrial, N)); x = np.log10(1 + zt)
    yy = yd(zt, zgrid, grid) + rng.normal(0, sigma_dex, size=(ntrial, N))
    xm = x - x.mean(1, keepdims=True); ym = yy - yy.mean(1, keepdims=True)
    s = (xm * ym).sum(1) / (xm * xm).sum(1); return s.mean(), s.std()

Ns = np.array([30, 50, 100, 200, 350, 500])
sepR, sepD = [], []
for N in Ns:
    sR, eR = forecast(N, gR); sD, eD = forecast(N, gD); sC, eC = forecast(N, gC)
    sepR.append(abs(sR - sC) / np.sqrt(eR ** 2 + eC ** 2))
    sepD.append(abs(sD - sC) / np.sqrt(eD ** 2 + eC ** 2))

fig, (axL, axR) = plt.subplots(1, 2, figsize=(11, 4.4))
zf = np.linspace(0, 3, 200)
axL.axhspan(-sigma_dex * 1e3, sigma_dex * 1e3, color="0.85", zorder=0, label="±0.05 dex per-object scatter")
axL.plot(zf, [1e3 * yR(z) for z in zf], color="C3", lw=2,
         label=r"rate $a_0\!\propto\! cH$ (+%.0f%% at $z{=}2$)" % (100 * (10 ** yR(2.0) - 1)))
d2 = 100 * (10 ** y(2.0, *DESI_REP) - 1)
for k in DESI:
    axL.plot(zf, [1e3 * y(z, *DESI[k]) for z in zf], color="C0", lw=1, alpha=0.5)
axL.plot(zf, [1e3 * y(z, *DESI_REP) for z in zf], color="C0", lw=2,
         label=r"density (DESI band, %.0f%% at $z_f{=}2$)" % d2)
axL.plot(zf, 0 * zf, color="k", lw=1.6, ls="--", label=r"constant $a_0$")
axL.set_xlabel(r"observation redshift $z_{\rm obs}$")
axL.set_ylabel(r"$\Delta\log_{10}\sigma$ at fixed $M_b$ (milli-dex)")
axL.set_xlim(0, 3); axL.set_ylim(-40, 190); axL.set_title("The dispersion observable")
axL.legend(fontsize=8.5, loc="upper left", framealpha=0.9)

axR.axhline(3.0, color="0.6", lw=1, ls=":"); axR.text(505, 3.2, r"$3\sigma$", color="0.4", fontsize=9, ha="right")
axR.plot(Ns, sepR, "o-", color="C3", lw=2, label="rate vs constant (killable)")
axR.plot(Ns, sepD, "s-", color="C0", lw=2, label="density vs constant")
axR.set_xlabel(r"sample size $N$ (isolated deep-MOND systems)")
axR.set_ylabel(r"separability from constant $a_0$ ($\sigma$)")
axR.set_xlim(20, 520); axR.set_ylim(0, 12); axR.set_title("Preregistered forecast")
axR.legend(fontsize=9, loc="upper left", framealpha=0.9)
fig.tight_layout()
save_fig(fig, "fig_forecast")
plt.close(fig)

# ---- print the exact numbers for the text --------------------------------
print("density a0(z=2)/a0(0) per DESI SN sample:")
for k in DESI:
    r = a0_density(2.0, *DESI[k]); print(f"  {k:10s} w0wa={DESI[k]}  a0(2)={r:.3f}  dsigma={100*(10**(0.25*np.log10(r))-1):+.1f}%")
print(f"SEDE example: a0(2)={a0_density(2.0,*SEDE):.3f}")
print(f"rate a0(2)={a0_rate(2.0):.3f}  dsigma={100*(10**yR(2.0)-1):+.1f}%")
print("separability vs N:  " + " ".join(f"N={n}:rate={r:.1f}σ,dens={d:.1f}σ" for n, r, d in zip(Ns, sepR, sepD)))

# ---- high-z (JWST epoch) a0 ratios and MOND growth boost -------------------
# Growth-RATE boost ~ a0^{1/4} (omega ~ sqrt(G_eff rho), G_eff ~ sqrt(a0)); the
# FORCE enhancement is sqrt(a0). MOND growth is nonlinear (Sanders 1998; McGaugh
# et al. 2024), so this linear estimate is illustrative; the ordering is robust.
print("\nJWST-epoch a0(z)/a0(0), growth boost B=a0^{1/4}, and PS abundance ratio [density = DES-Y5]:")
def abund_ratio(B, nu): return np.exp((nu**2/2)*(1 - 1/B**2))
for z in (8, 10):
    ad, ar = a0_density(z, *DESI_REP), a0_rate(z)
    Bd, Br = ad**0.25, ar**0.25
    print(f"  z={z}: a0 dens={ad:.2f} rate={ar:5.1f} | B dens={Bd:.2f} rate={Br:.2f}"
          f" | n/n_const @nu=4: dens={abund_ratio(Bd,4):.1e} rate={abund_ratio(Br,4):.1e}"
          f" @nu=5: dens={abund_ratio(Bd,5):.1e} rate={abund_ratio(Br,5):.1e}")
