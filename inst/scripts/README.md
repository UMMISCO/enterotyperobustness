# Scripts

None of these scripts is needed to use the package; they document how the
paper's numbers were, and will be, produced. They read data from paths given by
environment variables and write under `out/` (git-ignored).

| Script | Purpose | Needs |
|---|---|---|
| `reproduce_check.R` | Re-runs the stored `Healthy-Healthy` ERP experiment with the legacy seeds and compares draws, labels, maps, posteriors and flags with the 2022 files. | `ENTEROTYPES_MIRROR` |
| `metacardis_data.R` | `load_metacardis_data()` reads the master object (abundance, gene richness, clinical groups) the way the 2022 orchestrators did; `load_published_core_variable()` reads the Figure 1/2 core/variable lists. Sourced by the revision scripts. | `METACARDIS_DATA` (never published) |
| `run_metacardis.R` | *Planned (todo E6):* derives the 2022 anchors and pools and runs the 15 ERP conditions and the ERPM fractions with the package. | `METACARDIS_DATA` |
| `reproduce_figures.R` | Figures 1 to 4 and S2 to S4 from result objects. | outputs of the above or the mirror |
| `revision/01_…10_….R` | The reviewer-requested re-analyses, one per point; see `Reviews/todo.md` section A of the manuscript repository. | mirror (01 to 07) or master data (08 to 10) |
