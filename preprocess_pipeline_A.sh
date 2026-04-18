#!/bin/bash

set -euo pipefail

CHR="${1:-6}"
THREADS="${2:-4}"
SW_FORMULA="${3:-max_absz}"

WORKDIR="$HOME/nd_pipeline_2"
SCRIPTS="$WORKDIR/scripts"
DATA="$WORKDIR/data"
OUT="$WORKDIR/results/pipeline_A/$SW_FORMULA"

CHR_LEN_FILE="/home/share/human.data/ref.fa/hg19.chr.lengths/hg19.chrom.len"

NIS="$DATA/raw/IBS.YRI.grch37.chr${CHR}.em.tsv"
GTEX_DIR="$DATA/gtex/GTEx_Analysis_v7_eQTL"
GENCODE_GTF="$DATA/gencode/gencode.v19.annotation.gtf"
HAPMAP="$DATA/hapmap/genetic_map_GRCh37_chr${CHR}.txt"
CALLABLE_BED="$DATA/callability/chr${CHR}.renamed.bed"

mkdir -p "$OUT"

python3 "$SCRIPTS/01_prepare_windows_fw.py" \
  --chr-len-file "$CHR_LEN_FILE" \
  --chr "$CHR" \
  --nis "$NIS" \
  --outdir "$OUT"

python3 "$SCRIPTS/02_prepare_eqtl_sw.py" \
  --chr "$CHR" \
  --gtex-dir "$GTEX_DIR" \
  --windows-bed "$OUT/chr${CHR}_windows.bed" \
  --fw-tsv "$OUT/chr${CHR}_windows_Fw.tsv" \
  --outdir "$OUT" \
  --threads "$THREADS" \
  --sw-formula "$SW_FORMULA"

python3 "$SCRIPTS/03_prepare_tss_dtss.py" \
  --chr "$CHR" \
  --gtf "$GENCODE_GTF" \
  --windows-bed "$OUT/chr${CHR}_windows.bed" \
  --fw-sw-tsv "$OUT/chr${CHR}_windows_Fw_Sw.tsv" \
  --outdir "$OUT"

python3 "$SCRIPTS/04_prepare_recomb.py" \
  --chr "$CHR" \
  --hapmap "$HAPMAP" \
  --fw-sw-dtss-tsv "$OUT/chr${CHR}_windows_Fw_Sw_DTSS.tsv" \
  --outdir "$OUT"

python3 "$SCRIPTS/05_prepare_valid_windows.py" \
  --chr "$CHR" \
  --windows-bed "$OUT/chr${CHR}_windows.bed" \
  --callable-bed "$CALLABLE_BED" \
  --windows-full-tsv "$OUT/chr${CHR}_windows_full.tsv" \
  --outdir "$OUT" \
  --min-valid-frac 0.8
