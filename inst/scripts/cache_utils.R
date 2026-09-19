# Caches with a manifest. Every expensive intermediate result written by the
# revision scripts is stored as list(manifest, value); the manifest records the
# parameters (including an algorithm version string), an MD5 fingerprint of the
# actual inputs (the abundance matrix itself, sample sets, groups, or the 2022
# files' names, sizes and modification times), the package, DirichletMultinomial
# and R versions, the git commit and the date.
#
# A cache is reused only when the parameters, the input fingerprint and the
# DirichletMultinomial version match. A different package version or commit is
# reported, and refused when CACHE_STRICT is set (the archival run). A cache
# without manifest, or whose manifest differs, is recomputed, unless
# CACHE_MIGRATE is set: then the current manifest is written into it and it is
# used (a one-time local migration, never for an archival run).
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
cache_read <- function(file, manifest,
                       strict = nzchar(Sys.getenv("CACHE_STRICT")), migrate = nzchar(Sys.getenv("CACHE_MIGRATE"))) {
  if (!file.exists(file)) return(NULL)
  x <- readRDS(file)
  wrapped <- is.list(x) && !is.null(x$manifest) && "value" %in% names(x)
  adopt <- function(value, why) {
    if (!migrate) { message("cache ", file, ": ", why, "; recomputing (set CACHE_MIGRATE=1 to adopt it as is)"); return(NULL) }
    message("cache ", file, ": ", why, "; MIGRATED to the current manifest (CACHE_MIGRATE is set)")
    cache_write(value, file, manifest); value
  }
  if (!wrapped) return(adopt(x, "no manifest"))
  m <- x$manifest
  if (!identical(m$params, manifest$params)) return(adopt(x$value, "parameters differ"))
  if (!identical(m$inputs_md5, manifest$inputs_md5)) return(adopt(x$value, "input fingerprint differs"))
  if (!identical(m$dmm_version, manifest$dmm_version)) return(adopt(x$value, paste0("DirichletMultinomial ", m$dmm_version, " vs ", manifest$dmm_version)))
  other <- c(if (!identical(m$package_version, manifest$package_version)) paste0("package ", m$package_version, " vs ", manifest$package_version),
             if (!identical(m$git_commit, manifest$git_commit)) paste0("commit ", m$git_commit, " vs ", manifest$git_commit))
  if (length(other)) {
    if (strict) { message("cache ", file, ": ", paste(other, collapse = ", "), "; refused (CACHE_STRICT is set)"); return(NULL) }
    message("cache ", file, ": ", paste(other, collapse = ", "), " (accepted: parameters, inputs and DMM version match)")
  }
  message("cache ", file, " reused (written ", m$date, ", commit ", m$git_commit, ")")
  x$value
}
cache_write <- function(value, file, manifest) { saveRDS(list(manifest = manifest, value = value), file); invisible(value) }

# ---- manifests shared by several scripts -----------------------------------------------
# The 2022 ERPM environments read by script 00: fingerprint = names, sizes and mtimes of the files.
erpm_cache_manifest <- function(conf_dir) {
  files <- sort(list.files(conf_dir, pattern = "_TotalEnvironment\\.RData$", full.names = TRUE))
  fi <- file.info(files)
  cache_manifest(params = list(algorithm = "read_legacy_erpm v1: integer-coded assignments, stored Prob_Enterotype, flags, gap, ties"),
                 inputs = list(files = basename(files), size = fi$size, mtime = format(fi$mtime, "%Y-%m-%d %H:%M:%S")))
}
load_erpm_cache <- function(mirror, file = "out/cache/erpm_legacy.rds") {
  m <- erpm_cache_manifest(file.path(mirror, "analyses/2.Estimating_DMM/Article/Confidence"))
  x <- cache_read(file, m)
  if (is.null(x)) stop("No valid ERPM cache at ", file, ": run inst/scripts/revision/00_cache_erpm.R first.")
  x
}
# The full-cohort reference DMM fit of script 03 (also used by 08): fingerprint = the matrix itself.
reference_fit_manifest <- function(A, md) {
  cache_manifest(params = list(algorithm = "reference_centroids v1", K = 4, seed = 314, integerise = "trunc", annotation = "rank"),
                 inputs = list(A = A, marker_features = md$marker_features, richness = md$richness[rownames(A)]))
}
