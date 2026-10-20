# enterotyperobustness 1.0.0

First public release, archived at Zenodo with the revision of Creus-Martí *et al.*,
"A Sampling-Based Computational Framework for Evaluating Enterotype Assignment
Stability in Clinical Microbiome Studies" (*Scientific Reports*).

* ERP (composition perturbation) and ERPM (size perturbation) under DMM clustering,
  with rank and centroid-matching annotation of the enterotypes.
* Robustness, assignment entropy and core/variable classification; design-aware
  tests (Friedman/Conover, permutation chi-square, Cochran's Q, mixed models).
* Readers for the 2022 result files (`read_legacy_erp()`, `read_legacy_erpm()`).
* Scripts under `inst/scripts/` that produce the figures and tables of the paper
  from the private MetaCardis data (not distributed).
