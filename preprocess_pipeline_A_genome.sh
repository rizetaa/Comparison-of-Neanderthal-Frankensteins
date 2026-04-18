#!/bin/bash

set -euo pipefail

THREADS="${1:-4}"
SW_FORMULA="${2:-max_absz}"

WORKDIR="$HOME/nd_pipeline_2"
OUT="$WORKDIR/results/pipeline_A/$SW_FORMULA"

mkdir -p "$OUT"

for CHR in {1..22}; do
  bash scripts/preprocess_pipeline_A.sh "$CHR" "$THREADS" "$SW_FORMULA"
done

python3 - << PY
import glob
import pandas as pd
from pathlib import Path

workdir = Path("$WORKDIR")
out = Path("$OUT")
files = sorted(glob.glob(str(workdir / "data" / "raw" / "IBS.YRI.grch37.chr*.em.tsv")))
if not files:
    raise RuntimeError("No NIS files found")

lengths = []
for f in files:
    df = pd.read_csv(f, sep=r"\s+")
    if "Length" not in df.columns:
        raise RuntimeError(f"Column 'Length' not found in {f}")
    lengths.append(df["Length"])

all_lengths = pd.concat(lengths, ignore_index=True)
median_len = int(all_lengths.median())
(out / "median_seg_len.txt").write_text(str(median_len) + "\n", encoding="utf-8")
print(f"Genome-wide median NIS segment length: {median_len}")
print(f"Saved: {out / 'median_seg_len.txt'}")
PY

python3 scripts/merge_genomewide_windows.py --sw-formula "$SW_FORMULA"
