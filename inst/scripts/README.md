# Scripts

None of these scripts is needed to use the package; they document how the
paper's numbers were, and will be, produced. They read data from paths given by
environment variables and write under `out/` (git-ignored).

| Script | Purpose | Needs |
|---|---|---|
| `reproduce_check.R` | Re-runs the stored `Healthy-Healthy` ERP experiment with the legacy seeds and compares draws, labels, maps, posteriors and flags with the 2022 files. | `ENTEROTYPES_MIRROR` |
| `run_metacardis.R` | Derives anchors, pools, classes and gene richness from the MetaCardis master object exactly as the 2022 orchestrators did, and runs the 15 ERP conditions and the ERPM fractions. | `METACARDIS_DATA` (never published) |
| `reproduce_figures.R` | Figures 1 to 4 and S2 to S4 from result objects. | outputs of the above or the mirror |
| `revision/01_…10_….R` | The reviewer-requested re-analyses, one per point; see `Reviews/todo.md` section A of the manuscript repository. | mirror (01 to 07) or master data (08 to 10) |
