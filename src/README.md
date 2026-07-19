# src/ — APDM analysis code

The experiment scripts that reproduce every quantitative result in the paper
(`paper/APDM_galactic.md`). Run them from the repo root with `PYTHONPATH=src`, or
run all of them at once via `python3 reproduce.py` (from the repo root).

```bash
# single experiment
PYTHONPATH=src python3 src/branch_discrimination_test.py

# everything, with a PASS/FAIL summary
PYTHONPATH=src python3 reproduce.py
```

## Headline result (E-1): the +16.7σ density-branch exclusion

`branch_discrimination_test.py` is the headline test. It compares the two `a0(z)`
branches against the MUSE-DARK III rotator sample (Ciocan et al. 2026,
arXiv:2604.22613; `a0(z) = 1.00 + 1.59 z`, 0.33 < z < 1.44) using the DESI DR2
posterior templates in `data/desi_template_quantiles.csv`. It prints, branch-averaged:

* **DENSITY branch** `a0 ∝ √ρ_DE`: **+16.7σ — EXCLUDED** (falls with z; the data rise — wrong sign)
* **RATE branch** `a0 ∝ H(z)`: +4.4σ (right sign, amplitude short)

The +16.7σ is **conditional on the MUSE-DARK III modelling chain** (see §3/§8 of the paper).

## The `sede/` package ([SEDE] experiments)

The vendored `sede/` package (Structural-Entropy Dark Energy) supplies the
background/growth functions (`sede.friedmann`, `sede.theory`) that the cosmology
experiments import. It lives beside the scripts, so `PYTHONPATH=src` makes both the
scripts and `import sede` resolve. Only `numpy`/`scipy` are needed for the modules
used here; the package's optional heavy deps (camb, classy, emcee, …) are not required.

Experiments marked **[SEDE]** in `reproduce.py`:
`e9_mechanism_bridge.py`, `e9_tightened.py`, `solve_emergent.py`,
`compute_solidification.py`, `test_one_coupling_two_faces.py`,
`test_twofield_completion.py`, `test_doom_instability.py`,
`test_growth_perturbations.py`, `test_downstream_tasks.py`, `make_figures.py`.
Each also carries a graceful fallback if `sede` is unavailable.

## Network / [data] experiments

The following scripts read public catalogues. Each is **cache-guarded**: it uses the
data file already shipped in `data/` and only reaches the network **on first run if
that cache is missing**. Since the caches are vendored, they run offline as-is.

| script | catalogue | source |
|---|---|---|
| `test_law_vs_attractor_gema.py`  | GEMA-VAC 2.0.2 + Ristea+2024 | SDSS DR17 SAS / VizieR |
| `test_law_vs_attractor_manga.py` | MaNGA kinematics (Ristea+2024) | VizieR J/MNRAS/527/7438 |
| `test_reframe_rar_sparc.py`      | SPARC mass models (Lelli+2016) | astroweb.cwru.edu / GitHub mirror |
| `scaffold_law_vs_attractor.py`   | SPARC mass models | astroweb.cwru.edu / GitHub mirror |
| `test_wb_metallicity_iv.py`      | Gaia DR3 wide-binary self-join | ESA Gaia TAP |

`run_smoking_gun_predictions.py` and `run_screening_dwarf_test.py` describe
data-driven falsifiers (MaNGA/SPARC/SDSS); as shipped they are self-contained and do
not fetch, but a full observational run of those falsifiers requires network access
to the same catalogues on first use.

See `data/README.md` for the provenance and licensing of the shipped data files.
