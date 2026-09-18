# A small four-community count table with the marker structure the annotation
# expects: one Bacteroides-heavy low-diversity group (Bact2), one moderate
# Bacteroides group (Bact1), one Prevotella group, one Ruminococcus group.
simulate_counts <- function(n = 80, seed = 1, depth = 3000) {
  set.seed(seed)
  feats <- c("Bacteroides", "Prevotella", "Ruminococcus", "Faecalibacterium",
             "Subdoligranulum", "Eubacterium", paste0("g", 1:14))
  grp <- rep(1:4, length.out = n)
  base <- matrix(2, 4, length(feats), dimnames = list(NULL, feats))
  base[1, "Bacteroides"] <- 300                                   # Bact2
  base[2, "Bacteroides"] <- 120; base[2, "Faecalibacterium"] <- 60 # Bact1
  base[3, "Prevotella"] <- 300                                    # Prevotella
  base[4, "Ruminococcus"] <- 150; base[4, "Faecalibacterium"] <- 100
  base[4, "Subdoligranulum"] <- 80; base[4, "Eubacterium"] <- 40  # Ruminococcus
  X <- t(vapply(grp, function(g) {
    p <- stats::rgamma(length(feats), base[g, ]); p <- p / sum(p)
    stats::rmultinom(1, depth, p)[, 1] + stats::runif(length(feats), 0, 0.9)  # non-integer, like the real table
  }, numeric(length(feats))))
  dimnames(X) <- list(sprintf("s%03d", seq_len(n)), feats)
  attr(X, "group") <- grp
  attr(X, "diversity") <- stats::setNames(c(100, 200, 150, 300)[grp] + stats::rnorm(n, 0, 3), rownames(X))
  X
}
markers <- c(bacteroides = "Bacteroides", prevotella = "Prevotella")
