#!/bin/bash

set -euo pipefail

PLOT_TYPE="${1:-violin}"
SW_FORMULA="${2:-max_absz}"

WORKDIR="$HOME/nd_pipeline_2"
SCRIPTS="$WORKDIR/scripts"
PIPE_A="$WORKDIR/results/pipeline_A/$SW_FORMULA"
OUTBASE="$WORKDIR/results/analysis/genome/$SW_FORMULA"
QC_OUT="$OUTBASE/qc"
REG_OUT="$OUTBASE/regression"
STRAT_OUT="$OUTBASE/stratified_control"
ADAPT_OUT="$OUTBASE/adaptive_candidates"

mkdir -p "$QC_OUT" "$REG_OUT" "$STRAT_OUT" "$ADAPT_OUT"

python3 "$SCRIPTS/analysis_qc_introgression.py" \
  --input "$PIPE_A/genome_windows_full_valid.tsv" \
  --outdir "$QC_OUT" \
  --valid-only

python3 "$SCRIPTS/analysis_regression.py" \
  --input "$PIPE_A/genome_windows_full_valid.tsv" \
  --outdir "$REG_OUT" \
  --plot-type "$PLOT_TYPE"

python3 "$SCRIPTS/analysis_stratified_control.py" \
  --input "$PIPE_A/genome_windows_full_valid.tsv" \
  --outdir "$STRAT_OUT" \
  --plot-type "$PLOT_TYPE" \
  --median-seg-len-file "$PIPE_A/median_seg_len.txt" \
  --n-bootstrap 5000 \
  --threads 8

python3 "$SCRIPTS/merge_genomewide_lead_eqtl.py" \
  --sw-formula "$SW_FORMULA" \
  --workdir "$WORKDIR"

python3 "$SCRIPTS/analysis_adaptive_candidates.py" \
  --input "$PIPE_A/genome_windows_full_valid.tsv" \
  --lead-eqtl "$PIPE_A/genome_lead_eqtl.tsv" \
  --outdir "$ADAPT_OUT" \
  --fw-quantile 0.95 \
  --sw-quantile 0.95 \
  --max-gap-bp 20000
