#!/usr/bin/env python3

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Merge chr*_lead_eqtl.tsv into one genome-wide table."
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
        path = base / f"chr{chrom}_lead_eqtl.tsv"
        if not path.exists():
            missing.append(str(path))
            continue
            
        df = pd.read_csv(path, sep="\t", low_memory=False)
        # Нормализуем столбец
        if "chrom" not in df.columns:
            if "chr_win" in df.columns:
                df = df.rename(columns={"chr_win": "chrom"})
            elif "CHR" in df.columns:
                df = df.rename(columns={"CHR": "chrom"})
            elif "chr" in df.columns:
                df = df.rename(columns={"chr": "chrom"})
            else:
                df["chrom"] = str(chrom)
        df["chrom"] = df["chrom"].astype(str)
        tables.append(df)

    if missing:
        raise RuntimeError(
            "Some chromosome-level lead-eQTL files are missing:\n" + "\n".join(missing)
        )
        
    if not tables:
        raise RuntimeError(f"No chromosome lead-eQTL tables found in {base}")

    merged = pd.concat(tables, ignore_index=True)
    sort_cols = [c for c in ["chrom", "win_start", "win_end", "win_id"] if c in merged.columns]
    if sort_cols:
        merged = merged.sort_values(sort_cols, kind="mergesort").reset_index(drop=True)
        
    out_path = base / "genome_lead_eqtl.tsv"
    merged.to_csv(out_path, sep="\t", index=False)
    print(f"Merged rows: {len(merged)}")
    print(f"Saved: {out_path}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
