#!/bin/bash
# ------------------------------------------------------------------------------
# release_rerun.sh  Archival re-run of the revision analyses from a clean clone.
#
#   bash inst/scripts/release_rerun.sh <ENTEROTYPES_MIRROR> <METACARDIS_DATA> <N_CORES> <WORKDIR>
#
# Run from the package root of a working copy whose HEAD is committed. It:
#   1. clones that exact commit into <WORKDIR>/enterotyperobustness (must not exist);
#   2. checks the clone is clean, has no out/, and records commit and version;
#   3. records R, platform and package versions;
#   4. builds the package and runs R CMD check with the private data visible, so the
#      legacy tests run too;
#   5. runs inst/scripts/revision/run_all.sh with CACHE_STRICT=1;
#   6. writes <WORKDIR>/audit_report.txt.
#
# Everything stays on this computer. audit_report.txt holds only technical facts
# (commit, versions, status, durations, warning counts, error messages, file counts
# and checksums), no sample identifiers or results; read it before sharing it.
# ------------------------------------------------------------------------------
set -uo pipefail
[ $# -eq 4 ] || { sed -n 4,6p "$0"; exit 1; }
MIRROR="$1"; DATA="$2"; NCORES="$3"; WORK="$4"
[ -f DESCRIPTION ] && [ -d .git ] || { echo "run from the package root"; exit 1; }
[ -d "$MIRROR" ] || { echo "ENTEROTYPES_MIRROR is not a directory: $MIRROR"; exit 1; }
[ -f "$DATA" ]   || { echo "METACARDIS_DATA is not a file: $DATA"; exit 1; }
[ -z "$(git status --porcelain)" ] || { echo "working copy has uncommitted changes: commit first"; exit 1; }
SRC="$(pwd)"; COMMIT="$(git rev-parse HEAD)"
mkdir -p "$WORK"; WORK="$(cd "$WORK" && pwd)"; CLONE="$WORK/enterotyperobustness"
[ -e "$CLONE" ] && { echo "$CLONE already exists: choose an empty WORKDIR"; exit 1; }
R="$WORK/audit_report.txt"
md5f() { if command -v md5 >/dev/null; then md5 -q "$1"; else md5sum "$1" | cut -d' ' -f1; fi; }
say() { echo "$*" | tee -a "$R"; }

: > "$R"
say "# enterotyperobustness — archival re-run audit"
say "started:  $(date '+%Y-%m-%d %H:%M:%S %Z')"
say "machine:  $(uname -srm), $NCORES cores requested"
say ""
say "## 1. Clean clone"
git clone -q --no-local "$SRC" "$CLONE" && cd "$CLONE" && git checkout -q "$COMMIT" || { say "clone FAILED"; exit 1; }
say "commit:           $(git rev-parse HEAD)  ($(git log -1 --format='%cd' --date=short))"
say "subject:          $(git log -1 --format=%s)"
say "tags at commit:   $(git tag --points-at HEAD | tr '\n' ' ')"
say "version:          $(grep '^Version:' DESCRIPTION | cut -d' ' -f2)"
say "clean tree:       $([ -z "$(git status --porcelain --ignored)" ] && echo yes || echo NO)"
say "out/ absent:      $([ -e out ] && echo NO || echo yes)"
say "data file md5:    $(md5f "$DATA")  ($(basename "$DATA"))"
say ""
say "## 2. Environment"
Rscript -e 'cat("R:", R.version.string, "\n")
pk <- c("DirichletMultinomial","clue","testthat","pkgload","ggplot2","ggalluvial","patchwork","lme4","lmerTest","emmeans","mclogit","nnet","cluster","dplyr","tidyr","boot","dunn.test","PMCMRplus","knitr","rmarkdown")
for (p in pk) cat(sprintf("  %-22s %s\n", p, if (requireNamespace(p, quietly = TRUE)) as.character(packageVersion(p)) else "NOT INSTALLED"))
cat("  platform:", R.version$platform, "| BLAS:", tryCatch(extSoftVersion()[["BLAS"]], error = function(e) "?"), "\n")' 2>&1 | tee -a "$R"
say ""
say "## 3. R CMD build / check (private data visible: legacy tests run)"
cd "$WORK"; t0=$(date +%s)
R CMD build "$CLONE" > build.log 2>&1 && TGZ=$(ls enterotyperobustness_*.tar.gz | head -1) || { say "build FAILED (see $WORK/build.log)"; TGZ=""; }
if [ -n "$TGZ" ]; then
  ENTEROTYPES_MIRROR="$MIRROR" METACARDIS_DATA="$DATA" _R_CHECK_CRAN_INCOMING_REMOTE_=false R CMD check --no-manual "$TGZ" > check.log 2>&1
  say "tarball:          $TGZ (md5 $(md5f "$TGZ"))"
  say "check:            $(grep '^Status:' check.log)  [$(( $(date +%s) - t0 ))s]"
  grep -E '^\* checking .*(WARNING|ERROR|NOTE)$' check.log | sed 's/^/  /' | tee -a "$R"
  say "tests:            $(grep -hE '\[ FAIL' enterotyperobustness.Rcheck/tests/testthat.Rout* 2>/dev/null | tail -1)"
  grep -hA12 'Skipped tests' enterotyperobustness.Rcheck/tests/testthat.Rout* 2>/dev/null | grep -E '^\*' | sed 's/^/  skipped: /' | tee -a "$R"
fi
say ""
say "## 4. run_all.sh (CACHE_STRICT=1)"
cd "$CLONE"; t0=$(date +%s)
ENTEROTYPES_MIRROR="$MIRROR" METACARDIS_DATA="$DATA" CACHE_STRICT=1 N_CORES="$NCORES" bash inst/scripts/revision/run_all.sh > "$WORK/run_all.console.log" 2>&1
RC=$?
say "exit code:        $RC   total $(( ($(date +%s) - t0) / 60 )) min"
say ""
[ -f out/logs/run_all.log ] && sed 's/^/  /' out/logs/run_all.log | tee -a "$R"
say ""
say "per-script log: lines | warnings | errors | cache refused/recomputed"
for f in out/logs/*.log; do
  b=$(basename "$f" .log); [ "$b" = run_all ] && continue
  say "$(printf '  %-32s %6s | %4s | %4s | %4s' "$b" "$(wc -l < "$f" | tr -d ' ')" \
       "$(grep -ciE 'warning' "$f")" "$(grep -ciE '^error|error in|execution halted' "$f")" "$(grep -ciE 'refus|recomput|migrat' "$f")")"
done
if [ $RC -ne 0 ]; then
  failed=$(grep -o 'FAILED, see out/logs/[^ ]*' out/logs/run_all.log | tail -1 | sed 's/FAILED, see //')
  say ""; say "last 15 lines of $failed (REVIEW BEFORE SHARING: may quote data):"
  [ -n "$failed" ] && tail -15 "$failed" | sed 's/^/  | /' | tee -a "$R"
fi
say ""
say "## 5. Outputs"
if [ -f out/output_manifest.txt ]; then
  say "files in manifest: $(wc -l < out/output_manifest.txt | tr -d ' ')   manifest md5: $(md5f out/output_manifest.txt)"
  say "files per directory:"
  find out -type f ! -path 'out/logs/*' | awk -F/ '{print "  " $2}' | sort | uniq -c | tee -a "$R"
else
  say "no out/output_manifest.txt"
fi
say ""
say "finished: $(date '+%Y-%m-%d %H:%M:%S %Z')"
say "outputs:  $CLONE/out   (stay on this computer)"
