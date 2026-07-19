# src/data/ — third-party data (provenance & licensing)

These files are **public third-party catalogues**, cached here so the analysis
scripts reproduce offline. They are redistributed under their original terms and are
**not** covered by this repository's MIT/CC-BY licences. Cite the original sources if
you use them. Each is only re-fetched by its script if the cache is deleted.

| file | contents | source & citation | terms |
|---|---|---|---|
| `desi_template_quantiles.csv` | DESI DR2 posterior quantiles of `a0(z)/a0(0)` for the density and rate branches (per SN set: Pantheon+, DES-Y5, Union3), used by `branch_discrimination_test.py` | Derived in this program from DESI DR2 BAO + SN posteriors (DESI Collaboration 2025). Redistributed as a small derived summary table. | derived summary; public DESI data |
| `GEMA_2.0.2.fits` | GEMA-VAC 2.0.2 galaxy-environment value-added catalogue (tidal strength Q_nn, close-pair merger probability) | Argudo-Fernández et al., SDSS DR17 GEMA-VAC. Original: `https://data.sdss.org/sas/dr17/manga/gema/2.0.2/GEMA_2.0.2.fits` | SDSS data policy |
| `MaNGA_Ristea2024.tsv`, `MaNGA_Ristea2024_id.tsv` | MaNGA DR17 stellar/gas kinematics (V/σ, TFR residuals) | Ristea et al. 2024, MNRAS 527, 7438. VizieR `J/MNRAS/527/7438` | VizieR / CDS terms |
| `SPARC_MassModels.mrt` | SPARC rotation-curve mass models | Lelli, McGaugh & Schombert 2016, AJ 152, 157. Official: `http://astroweb.cwru.edu/SPARC/` (mirror: `github.com/jobovy/sparc-rotation-curves`) | SPARC terms of use |
| `wb_gaia_pairs.csv` | Gaia DR3 wide-binary pairs with GSP-Phot [M/H] (d < 50 pc self-join), for the metallicity-IV test | ESA Gaia DR3 (Gaia Collaboration 2023), queried via the Gaia TAP service | ESA Gaia data licence (CC-BY-like) |

**How to regenerate:** delete the relevant cache file and re-run its script with
network access; the script will re-download from the source URL above.
