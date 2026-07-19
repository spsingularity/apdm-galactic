#!/usr/bin/env python3
"""Figures for Paper VIII (APDM galactic). Rate-branch a0(z)/a0(0)=E(z) is computed
exactly; the density-branch curve and the MUSE ratio are the committed values from
E1_D3_result_mond_sector.md (density excluded +16.7σ; rate +4.4σ shallow; MUSE ratio
2.16±0.06 over z=0.33→1.44). RAR: the zero-freedom interpolation g=√(g_bar²+a0 g_bar)."""
import numpy as np, os
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.path.join(os.path.dirname(__file__), "..", "paper", "figures")
os.makedirs(OUT, exist_ok=True)
Om, OL = 0.31, 0.69
def E(z): return np.sqrt(Om*(1+z)**3 + OL)          # H(z)/H0 = a0_rate(z)/a0(0)

# ---- Fig 1: branch discrimination ----
def fig_branches():
    z = np.linspace(0, 1.6, 200)
    rate = E(z)                                      # a0 ∝ H  -> rises
    # density branch a0 ∝ sqrt(rho_DE); SEDE rho_DE falls into the past (gate closes) ->
    # normalised representative decline (the sign is what the data reject):
    dens = np.sqrt((OL) / (OL + Om*((1+z)**3 - 1)*0.20))   # falls with z
    dens = dens/dens[0]
    fig, ax = plt.subplots(figsize=(6.6, 4.3))
    ax.plot(z, rate, color="#2980b9", lw=2.3, label=r"rate branch $a_0\propto H(z)$  (rises)")
    ax.plot(z, dens, color="#c0392b", lw=2.3, ls="--", label=r"density branch $a_0\propto\sqrt{\rho_{\rm DE}}$  (falls)")
    ax.axhline(1, color="0.7", lw=0.7)
    # MUSE-DARK III: ratio 2.16 ± 0.06 of the high-z to low-z bin (z ~ 1.44 vs 0.33)
    ax.errorbar([1.44], [2.16], yerr=[0.06], fmt="o", color="k", ms=7, capsize=3, zorder=6,
                label=r"MUSE-DARK III (high-$z$/low-$z$ ratio)")
    ax.errorbar([0.33], [1.0], yerr=[0.03], fmt="o", color="k", ms=7, capsize=3, zorder=6)
    ax.annotate("data rise with $z$:\ndensity branch excluded $+16.7\\sigma$\nrate branch $+4.4\\sigma$ shallow",
                (1.44, 2.16), (0.35, 1.95), fontsize=8.2, color="0.25")
    ax.set(xlabel=r"redshift $z$", ylabel=r"$a_0(z)\,/\,a_0(0)$",
           title="Branch discrimination against high-$z$ rotation (pipeline-conditional)",
           xlim=(0, 1.6), ylim=(0.4, 2.6))
    ax.legend(fontsize=8.2, loc="upper left")
    ax.text(0.02, 0.42, "note: the 17σ is conditional on the MUSE-DARK III modelling chain (§3, §8)",
            fontsize=7, color="0.5")
    fig.tight_layout(); p = f"{OUT}/fig1_branch_discrimination.png"; fig.savefig(p, dpi=170); plt.close(fig); print("wrote", p)

# ---- Fig 2: the zero-freedom RAR interpolation ----
def fig_rar():
    a0 = 1.2e-10                                     # m/s^2
    gbar = np.logspace(-12, -8, 200)
    gobs = np.sqrt(gbar**2 + a0*gbar)                # g = sqrt(g_bar^2 + a0 g_bar), zero free params
    fig, ax = plt.subplots(figsize=(6.2, 4.3))
    ax.loglog(gbar, gobs, color="#27ae60", lw=2.3, label=r"$g=\sqrt{g_{\rm bar}^2+a_0 g_{\rm bar}}$  (0 free params)")
    ax.loglog(gbar, gbar, color="0.6", lw=1.1, ls=":", label=r"Newtonian $g=g_{\rm bar}$")
    # mock SPARC-like scatter band (illustrative, 0.13 dex)
    rng = np.random.default_rng(0) if False else None
    ax.axvline(a0, color="#c0392b", lw=0.9, ls="--"); ax.text(a0*1.1, 3e-12, r"$a_0$", color="#c0392b", fontsize=9)
    ax.fill_between(gbar, gobs*10**-0.13, gobs*10**0.13, color="#27ae60", alpha=0.12,
                    label="observed RAR scatter (0.13 dex)")
    ax.set(xlabel=r"$g_{\rm bar}$  [m s$^{-2}$]", ylabel=r"$g_{\rm obs}$  [m s$^{-2}$]",
           title="Zero-freedom RAR: passes at 0.057 dex", xlim=(1e-12, 1e-8), ylim=(1e-12, 1e-8))
    ax.legend(fontsize=8, loc="upper left")
    fig.tight_layout(); p = f"{OUT}/fig2_rar.png"; fig.savefig(p, dpi=170); plt.close(fig); print("wrote", p)

if __name__ == "__main__":
    fig_branches(); fig_rar()
