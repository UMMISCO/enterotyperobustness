#' Fit one Dirichlet Multinomial Mixture model
#'
#' Thin wrapper around [DirichletMultinomial::dmn()] that makes the two choices
#' the 2022 analyses left implicit explicit: how non-integer cells are turned
#' into counts, and the internal seed of the EM initialisation.
#'
#' @param counts Numeric matrix, samples in rows, features in columns.
#' @param K Number of mixture components.
#' @param seed Seed passed to `dmn()`. The published analyses used 314 for
#'   every iteration.
#' @param integerise How to turn non-integer cells into counts; see
#'   [integerise_counts()]. `"trunc"` reproduces the published fits.
#' @param verbose Passed to `dmn()`.
#' @return A `DirichletMultinomial::DMN` object.
#' @seealso [dmm_alpha()], [dmm_posterior()], [dmm_assign()]
#' @importFrom DirichletMultinomial dmn mixture fitted
#' @importFrom stats setNames sd
#' @examples
#' X <- simulate_enterotype_counts(n = 40)
#' fit <- fit_dmm(X, K = 4, seed = 314)
#' table(dmm_assign(fit))
#' @export
fit_dmm <- function(counts, K = 4, seed = 314,
                    integerise = c("trunc", "round", "none"), verbose = FALSE) {
  counts <- check_abundance(counts)
  counts <- integerise_counts(counts, match.arg(integerise))
  DirichletMultinomial::dmn(counts, k = K, seed = seed, verbose = verbose)
}

#' Dirichlet parameters of a fit, features x components
#' @param fit A `DMN` object.
#' @export
dmm_alpha <- function(fit) DirichletMultinomial::fitted(fit)

#' Component centroids of a fit: alpha normalised to sum to one per component
#' @param fit A `DMN` object.
#' @return A features x components matrix of expected proportions.
#' @export
dmm_centroids <- function(fit) {
  a <- dmm_alpha(fit)
  sweep(a, 2, colSums(a), "/")
}

#' Posterior component probabilities of a fit, samples x components
#' @param fit A `DMN` object.
#' @export
dmm_posterior <- function(fit) DirichletMultinomial::mixture(fit, assign = FALSE)

#' Maximum a posteriori component of every sample
#' @param fit A `DMN` object.
#' @return An integer vector, one component index per sample.
#' @export
dmm_assign <- function(fit) as.integer(DirichletMultinomial::mixture(fit, assign = TRUE))
