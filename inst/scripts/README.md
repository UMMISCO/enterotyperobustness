# Scripts

None of these scripts is needed to use the package; they document how the
paper's numbers were, and will be, produced. They read data from paths given by
environment variables (`ENTEROTYPES_MIRROR`: the archive of the 2022 analyses;
`METACARDIS_DATA`: the master data object; both private, neither has a default)
and write under `out/` (git-ignored). Run them from the package root
(`Rscript inst/scripts/revision/01_paired_tests_ERP.R`): they call
`pkgload::load_all(".")` and use relative `out/` paths.

| Script | Purpose | Needs |
|---|---|---|
| `reproduce_check.R` | Re-runs the stored `Healthy-Healthy` ERP experiment with the legacy seeds and compares draws, labels, maps, posteriors and flags with the 2022 files. | `ENTEROTYPES_MIRROR` |
| `metacardis_data.R` | `load_metacardis_data()` reads the master object (abundance, gene richness, clinical groups) the way the 2022 orchestrators did; `load_published_core_variable()` reads the Figure 1/2 core/variable lists. Sourced by the revision scripts. | `METACARDIS_DATA` (never published) |
| `run_metacardis.R` | *Planned (todo E6):* derives the 2022 anchors and pools and runs the 15 ERP conditions and the ERPM fractions with the package. | `METACARDIS_DATA` |
| `reproduce_figures.R` | *Planned (todo E7):* Figures 1 to 4 and S2 to S4 from result objects (panels B to F; the flower panels A and D are added in Inkscape). | outputs of the above or the mirror |
| `determinism_check.R`, `determinism_full.R`, `rerun_compare.R` | Run-to-run (non-)determinism of `dmn()` on the stored draws; see the note in `erp()`. | `ENTEROTYPES_MIRROR` |

## Revision scripts (`revision/`)

One per reviewer point; see `Reviews/todo.md` section A of the manuscript
repository. Run them in this order, as the later ones read the outputs of the
earlier ones under `out/`:

| Script | Purpose | Needs |
|---|---|---|
| `00_cache_erpm.R` | Reads the 56 stored ERPM environments once into `out/cache/erpm_legacy.rds`. | mirror |
| `01_paired_tests_ERP.R` | Friedman/Conover, permutation Chi-square and Cochran's Q on the four ERP designs (R1 major 5). | mirror |
| `01b_mixed_models_ERP.R` | Linear mixed model of entropy and mixed multinomial logit of labels; regenerated panels 1C/1F/2C/2F and Table S4. | mirror, `out/01_paired_tests` |
| `02_bootstrap_CI_ERPM.R` | BCa bootstrap ribbons for Figures 3, 4 and S3, next to the submitted and the corrected normal intervals (R1 major 6). | `out/cache`, mirror figure objects |
| `03_cluster_matching.R` | Scripts 03 and 04: centroid matching against a reference fit, adjusted Rand index, Bacteroides gap (R1 major 1). | mirror, `METACARDIS_DATA` |
| `05_anomaly_tables.R` | Tables S1 to S3 from the raw flags; re-aggregation without flagged iterations. Set `SUPPLEMENTAL_TEX` to compare with the submitted tables. | mirror, `out/cache` |
| `06_core_variable_sensitivity.R` | Unanimity vs all-but-one vs SD thresholds for the core/variable rule (R1 minor 1). | mirror, `out/cache` |
| `07` to `10` | *Planned:* ERPM baseline at the ERP sizes, ERP rerun with matching, PAM on Jensen-Shannon, alternative K. | mirror or `METACARDIS_DATA` |
