#!/bin/bash
# Runs the revision analyses 00 to 10 in order, from the package root, with one
# log per script under out/logs/ and an MD5 manifest of every output file.
#
#   ENTEROTYPES_MIRROR=/path/to/ENTEROTYPES METACARDIS_DATA=/path/to/mcdb.sub.26122022.Rda \
#   CACHE_STRICT=1 N_CORES=12 bash inst/scripts/revision/run_all.sh
#
# CACHE_STRICT=1 refuses any cache written by another commit or package version
# (use it for the archival run, from a clean clone without out/). SUPPLEMENTAL_TEX
# is optional (script 05 compares with the submitted tables when it is set).
set -euo pipefail
cd "$(dirname "$0")/../../.."
[ -f DESCRIPTION ] || { echo "run from the package root"; exit 1; }
: "${ENTEROTYPES_MIRROR:?set ENTEROTYPES_MIRROR (private archive of the 2022 analyses)}"
: "${METACARDIS_DATA:?set METACARDIS_DATA (private master data object)}"
[ -d "$ENTEROTYPES_MIRROR" ] || { echo "ENTEROTYPES_MIRROR is not a directory"; exit 1; }
[ -f "$METACARDIS_DATA" ] || { echo "METACARDIS_DATA is not a file"; exit 1; }
N_CORES="${N_CORES:-8}"
mkdir -p out/logs
scripts=(00_cache_erpm 01_paired_tests_ERP 01b_mixed_models_ERP 02_bootstrap_CI_ERPM 03_cluster_matching
         05_anomaly_tables 06_core_variable_sensitivity 07_size_baseline 08_erp_rerun_matched
         09_alt_clustering_PAM_JSD 10_alt_K)
echo "run_all started $(date) on commit $(git rev-parse --short HEAD 2>/dev/null || echo unknown)" | tee out/logs/run_all.log
for s in "${scripts[@]}"; do
  echo "== $s  $(date)" | tee -a out/logs/run_all.log
  if Rscript "inst/scripts/revision/$s.R" "$N_CORES" > "out/logs/$s.log" 2>&1; then
    echo "   ok" | tee -a out/logs/run_all.log
  else
    echo "   FAILED, see out/logs/$s.log" | tee -a out/logs/run_all.log; exit 1
  fi
done
find out -type f ! -path "out/logs/*" ! -name output_manifest.txt -print0 | sort -z | xargs -0 md5 -r 2>/dev/null > out/output_manifest.txt \
  || find out -type f ! -path "out/logs/*" ! -name output_manifest.txt -print0 | sort -z | xargs -0 md5sum > out/output_manifest.txt
echo "run_all finished $(date); $(wc -l < out/output_manifest.txt) output files listed in out/output_manifest.txt" | tee -a out/logs/run_all.log
