#!/bin/bash

set -euo pipefail

CHR="${1:-6}"
PLOT_TYPE="${2:-violin}"
SW_FORMULA="${3:-max_absz}"

WORKDIR="$HOME/nd_pipeline"
SCRIPTS="$WORKDIR/scripts"
PIPE_A="$WORKDIR/results/pipeline_A/$SW_FORMULA"
OUTBASE="$WORKDIR/results/analysis/chr${CHR}/$SW_FORMULA"

QC_OUT="$OUTBASE/qc"
REG_OUT="$OUTBASE/regression"
STRAT_OUT="$OUTBASE/stratified_control"
ADAPT_OUT="$OUTBASE/adaptive_candidates"

mkdir -p "$QC_OUT" "$REG_OUT" "$STRAT_OUT" "$ADAPT_OUT"

#python3 "$SCRIPTS/analysis_qc_introgression.py" \
#  --input "$PIPE_A/chr${CHR}_windows_full_valid.tsv" \
#  --outdir "$QC_OUT" \
#  --valid-only

#python3 "$SCRIPTS/analysis_regression.py" \
#  --input "$PIPE_A/chr${CHR}_windows_full_valid.tsv" \
#  --outdir "$REG_OUT" \
#  --plot-type "$PLOT_TYPE"

#python3 "$SCRIPTS/analysis_stratified_control.py" \
#  --input "$PIPE_A/chr${CHR}_windows_full_valid.tsv" \
#  --outdir "$STRAT_OUT" \
#  --plot-type "$PLOT_TYPE" \
#  --median-seg-len-file "$PIPE_A/chr${CHR}_median_seg_len.txt" \
#  --n-bootstrap 10000 \
#  --threads 8

python3 "$SCRIPTS/analysis_adaptive_candidates.py" \
  --input "$PIPE_A/chr${CHR}_windows_full_valid.tsv" \
  --lead-eqtl "$PIPE_A/chr${CHR}_lead_eqtl.tsv" \
  --outdir "$ADAPT_OUT" \
  --fw-quantile 0.95 \
  --sw-quantile 0.95 \
  --max-gap-bp 0
