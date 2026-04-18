#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge chr*_windows_full.tsv into one genome-wide table."
    )
    parser.add_argument(
        "--sw-formula",
        required=True,
        help="Sw formula subdirectory, e.g. max_absz, mean_absz, sum_absz",
    )
    parser.add_argument(
        "--workdir",
        default="~/nd_pipeline_2",
        help="Project root directory (default: ~/nd_pipeline_2)",
    )
    return parser.parse_args()

def main() -> int:
    args = parse_args()
    workdir = Path(args.workdir).expanduser()
    base = workdir / "results" / "pipeline_A" / args.sw_formula
    if not base.exists():
        raise RuntimeError(f"Directory does not exist: {base}")

    tables = []
    missing = []
    for chrom in range(1, 23):
        path = base / f"chr{chrom}_windows_full_valid.tsv"
        if not path.exists():
            missing.append(str(path))
            continue

        df = pd.read_csv(path, sep="\t")
        if "chrom" not in df.columns:
            df["chrom"] = str(chrom)
        else:
            df["chrom"] = df["chrom"].astype(str)
        tables.append(df)
        
    if missing:
        raise RuntimeError(
            "Some chromosome-level files are missing:\n" + "\n".join(missing)
        )
    if not tables:
        raise RuntimeError(f"No chromosome tables found in {base}")

    merged = pd.concat(tables, ignore_index=True)
    # Сортировка
    sort_cols = [c for c in ["chrom", "win_start", "win_end", "win_id"] if c in merged.columns]
    if sort_cols:
        merged = merged.sort_values(sort_cols, kind="mergesort").reset_index(drop=True)
    out_path = base / "genome_windows_full_valid.tsv"
    merged.to_csv(out_path, sep="\t", index=False)
    print(f"Merged rows: {len(merged)}")
    print(f"Saved: {out_path}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
