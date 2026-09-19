# Caches with a manifest. Every expensive intermediate result written by the
# revision scripts is stored as list(manifest, value); the manifest records the
# parameters, a fingerprint of the inputs, the package, DirichletMultinomial and
# R versions, the git commit and the date. A cache is reused only when the
# parameters and the input fingerprint match; otherwise it is recomputed.
# Caches written before this convention (19 Sep 2026, raw objects without a
# manifest) are adopted with a warning so that work in progress is not lost;
# the archived deposit must be produced from a clean clone without caches.
cache_fingerprint <- function(data) {
  tf <- tempfile(); on.exit(unlink(tf)); saveRDS(data, tf, version = 3); unname(tools::md5sum(tf))
}
cache_manifest <- function(params, inputs = NULL) {
  pv <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
  list(params = params, inputs_md5 = if (is.null(inputs)) NA_character_ else cache_fingerprint(inputs),
       package_version = pv("enterotyperobustness"), dmm_version = pv("DirichletMultinomial"),
       r_version = R.version.string,
       git_commit = tryCatch(suppressWarnings(system2("git", c("rev-parse", "--short", "HEAD"), stdout = TRUE, stderr = FALSE))[1], error = function(e) NA_character_),
       date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"))
}
cache_read <- function(file, manifest, adopt = TRUE) {
  if (!file.exists(file)) return(NULL)
  x <- readRDS(file)
  if (!(is.list(x) && !is.null(x$manifest) && !is.null(x$value))) {
    if (!adopt) { message("cache ", file, " has no manifest: ignored"); return(NULL) }
    message("cache ", file, " has no manifest (written before the manifest convention): adopted as is; ",
            "recompute from a clean clone before archiving")
    return(x)
  }
  m <- x$manifest
  if (!identical(m$params, manifest$params)) { message("cache ", file, ": parameters differ, recomputing"); return(NULL) }
  if (!identical(m$inputs_md5, manifest$inputs_md5)) { message("cache ", file, ": inputs differ, recomputing"); return(NULL) }
  if (!identical(m$dmm_version, manifest$dmm_version)) message("cache ", file, ": written with DirichletMultinomial ", m$dmm_version, ", now ", manifest$dmm_version)
  message("cache ", file, " reused (written ", m$date, ", commit ", m$git_commit, ")")
  x$value
}
cache_write <- function(value, file, manifest) { saveRDS(list(manifest = manifest, value = value), file); invisible(value) }
