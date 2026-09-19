# enterotyperobustness

Sampling-based evaluation of enterotype assignment stability under Dirichlet
Multinomial Mixture (DMM) clustering: the **ERP** (composition perturbation) and
**ERPM** (size perturbation) procedures of Creus-Martí *et al.* (*Scientific
Reports*, in revision). The package is data-agnostic: it takes an abundance
matrix and sample-ID vectors, and ships no cohort data.

```r
library(enterotyperobustness)
markers <- c(bacteroides = "Bacteroides", prevotella = "Prevotella")

# ERP: hold `anchors` fixed, add `n_complementary` samples from `pool` each iteration
res <- erp(abundance, anchors, pool, n_complementary = 86, n_iter = 500,
           marker_features = markers, diversity = gene_richness, n_cores = 10)
freq <- assignment_frequency(res)     # samples x 4, over the iterations drawn
rob  <- robustness(res)               # consensus enterotype and its frequency
ent  <- assignment_entropy(freq)

# ERPM: draw a stratified fraction of the cohort each iteration
res_m <- erpm(abundance, classes = clinical_group, n_iter = 500,
              marker_features = markers, diversity = gene_richness, n_cores = 10)

# Core / variable across ERP conditions (published rule: unanimity)
lab <- labels_by_condition(list(HH = res_hh, HO = res_ho, HD = res_hd))
cv  <- core_variable(lab, rule = "unanimous")
```

## Reproducing the published analyses

The defaults reproduce the 2022 runs: anchor drawn with `set.seed(314)` and 50 %
of its cohort (done by the caller; the MetaCardis loader is
`inst/scripts/metacardis_data.R`, and the full re-run script is planned), one
draw per iteration with `set.seed(i)`, `dmn(..., seed = 314)`, non-integer cells
**truncated** (`integerise = "trunc"`, what DirichletMultinomial 1.28 did
silently), rank annotation. `n_complementary` is the integer the original
computed as `round(0.3 * cohort_size)`: 86, 135 and 125 for the healthy,
T2D-metformin and obese anchors.

Two annotation methods: `annotate = "rank"` (published: Bact2 = highest mean
*Bacteroides*, Bact1 = second, Prevotella = highest *Prevotella*, Ruminococcus
= the rest) and `annotate = "match"` (Hungarian matching of component centroids
to a reference by Jensen-Shannon distance; needs `reference = reference_centroids(fit)`
of a reference fit). Annotation anomalies are returned as flags, never as
errors; the Bact1/Bact2 gap is recorded at every iteration.

`simulate_enterotype_counts()` gives a toy table so that every example, test and
the vignette (`vignette("enterotyperobustness")`) run without cohort data.
Seeded functions restore the caller's random stream on exit.

`read_legacy_erp()` and `read_legacy_erpm()` convert the 2022 result files into
package objects, so the reviewer-requested re-analyses run on the archived fits
without refitting.

## Installation

```r
# install.packages("BiocManager"); BiocManager::install("DirichletMultinomial")
devtools::install()     # from the package directory
```

## Repositories

- Source of truth: `https://git.ummisco.fr/jdzucker/enterotyperobustness` (UMMISCO GitLab, private during the revision).
- Mirror: `https://github.com/UMMISCO/enterotyperobustness`. Until GitLab push mirroring is confirmed to work, `git push` from the package pushes to both (the `origin` remote has two push URLs).
- No data is ever committed: the MetaCardis tables are read from paths given by environment variables (`METACARDIS_DATA`, `ENTEROTYPES_MIRROR`), which have no default. Without them the scripts under `inst/scripts/` stop with a message and the legacy tests are skipped.
