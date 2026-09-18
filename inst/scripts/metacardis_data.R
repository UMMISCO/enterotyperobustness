# Load the MetaCardis objects the way the 2022 orchestrators did. Not part of
# any published result: the path comes from METACARDIS_DATA and the file never
# enters the repository.
#
#   d <- load_metacardis_data()
#   d$abundance   2 022 samples x 117 genera (mOTU, rarefied to 10 M reads, "scaled down10",
#                 non-integer; the "-1" unclassified row dropped as in 2022)
#   d$richness    gene richness at 10 M reads, named by sample
#   d$groups      clinical group PATGROUPFINAL_C (column 162 of the bioclinical table), named by sample
#   d$bioclin     the full bioclinical table (2 214 rows; 2 022 have mOTU data)
load_metacardis_data <- function(path = Sys.getenv("METACARDIS_DATA",
                                                   "/Users/jdz/PROGRAMMATION/R/ENTEROTYPES/data/mcdb.sub.26122022.Rda")) {
  e <- new.env(); load(path, envir = e)
  m <- e$mcdb.sub
  g <- m$metagenomics$taxcomp$motu$scaled$down10$genus
  abundance <- t(as.matrix(g[rownames(g) != "-1", ]))
  bioclin <- m$bioclin$all.bioclin
  groups <- stats::setNames(as.character(bioclin[["PATGROUPFINAL_C"]]), rownames(bioclin))   # column 162 in 2022
  groups <- groups[names(groups) %in% rownames(abundance)]
  richness <- stats::setNames(m$metagenomics$comm$genecount[, "down_10000000"], rownames(m$metagenomics$comm$genecount))
  missing <- setdiff(rownames(abundance), names(richness))
  if (length(missing)) stop(length(missing), " sample(s) have no gene richness, e.g. ", missing[1])
  list(abundance = abundance, richness = richness[rownames(abundance)], groups = groups, bioclin = bioclin,
       marker_features = c(bacteroides = "Bacteroides", prevotella = "Prevotella"),
       extra_checks = list(FS = c("Faecalibacterium", "Subdoligranulum"),
                           FSR = c("Faecalibacterium", "Subdoligranulum", "Ruminococcus"),
                           FSRE = c("Faecalibacterium", "Subdoligranulum", "Ruminococcus", "Eubacterium")))
}

# The published core / variable lists (Figures 1 and 2 scripts, unanimity rule),
# one entry per (anchor set, sample, enterotype). `fig` is the
# 9.analyses_mars2026 directory of the mirror.
load_published_core_variable <- function(fig) {
  get_list <- function(rel, slot) { e <- new.env(); load(file.path(fig, rel), envir = e); e$samplesList_alluvial[[slot]] }
  sl <- list(healthy = get_list("Figure1/samplesList_alluvials_healthyCore.rda", "five_experiments"),
             t2d = get_list("Figure2/Figure2ABC_T2Dcore/samplesList_alluvials_T2DmetfCore.rda", "four_experiments"),
             obese = get_list("Figure2/Figure2DEF_ObeseCore/samplesList_alluvials_ObeseCore.rda", "four_experiments"))
  tab <- do.call(rbind, lapply(names(sl), function(a) do.call(rbind, lapply(names(sl[[a]]), function(e) rbind(
    data.frame(anchor_set = a, sample = sl[[a]][[e]]$core, enterotype = e, status = "core", stringsAsFactors = FALSE),
    data.frame(anchor_set = a, sample = sl[[a]][[e]]$variable, enterotype = e, status = "variable", stringsAsFactors = FALSE))))))
  list(lists = sl, table = tab)
}
