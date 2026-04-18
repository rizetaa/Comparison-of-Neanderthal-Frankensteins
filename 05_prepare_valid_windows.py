#!/usr/bin/env python3

from __future__ import annotations
from pathlib import Path

import argparse
import shutil
import subprocess
import sys
import pandas as pd

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Annotate windows with callable/valid fraction and merge into windows_full table."
    )
    parser.add_argument(
        "--chr",
        required=True,
        dest="chrom",
        help="Chromosome name, e.g. 6 or X",
    )
    parser.add_argument(
        "--windows-bed",
        required=True,
        help="Path to chr*_windows.bed",
    )
    parser.add_argument(
        "--callable-bed",
        required=True,
        help="BED with callable/valid regions for this chromosome (BED3/BED4/BEDn)",
    )
    parser.add_argument(
        "--windows-full-tsv",
        required=True,
        help="Path to chr*_windows_full.tsv",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory",
    )
    parser.add_argument(
        "--min-valid-frac",
        type=float,
        default=0.8,
        help="Minimum callable fraction to mark window as valid (default: 0.8)",
    )
    return parser.parse_args()

def check_dependencies() -> None:
    if shutil.which("bedtools") is None:
        raise RuntimeError("bedtools not found in PATH")

def load_windows(windows_bed: Path, chrom: str) -> pd.DataFrame:
    df = pd.read_csv(
        windows_bed,
        sep="\t",
        header=None,
        names=["chrom", "win_start", "win_end", "win_id"],
        dtype={"chrom": str, "win_start": int, "win_end": int, "win_id": int},
    )
    df = df[df["chrom"].astype(str) == str(chrom)].copy()
    if df.empty:
        raise RuntimeError(f"No windows found for chromosome {chrom} in {windows_bed}")

    df["win_len"] = df["win_end"] - df["win_start"]
    return df

def run_intersect(windows_bed: Path, callable_bed: Path, out_path: Path) -> None:
    cmd = [
        "bedtools",
        "intersect",
        "-a",
        str(windows_bed),
        "-b",
        str(callable_bed),
        "-wo",
    ]
    with out_path.open("w", encoding="utf-8") as out:
        subprocess.run(cmd, stdout=out, check=True)

def compute_validity(
    windows_df: pd.DataFrame,
    intersect_path: Path,
    min_valid_frac: float,
) -> pd.DataFrame:
    if intersect_path.stat().st_size == 0:
        validity = windows_df[["win_id", "chrom", "win_start", "win_end", "win_len"]].copy()
        validity["valid_bp"] = 0
        validity["valid_frac"] = 0.0
        validity["is_valid"] = 0
        return validity

    ix = pd.read_csv(intersect_path, sep="\t", header=None, low_memory=False)
    # A = окна BED: chrom, win_start, win_end, win_id
    # B = callable BED: мб >=3 c столбцов
    # последний столбец = overlap_bp
    if ix.shape[1] < 8:
        raise RuntimeError(
            f"Unexpected intersect format in {intersect_path}: got {ix.shape[1]} columns"
        )
      
    ix = ix.rename(
        columns={
            0: "chrom_win",
            1: "win_start",
            2: "win_end",
            3: "win_id",
            ix.shape[1] - 1: "overlap_bp",
        }
    )
    ix["win_id"] = pd.to_numeric(ix["win_id"], errors="raise").astype(int)
    ix["overlap_bp"] = pd.to_numeric(ix["overlap_bp"], errors="raise").astype(int)
    valid_bp = (
        ix.groupby("win_id", as_index=False)["overlap_bp"]
        .sum()
        .rename(columns={"overlap_bp": "valid_bp"})
    )
    validity = windows_df[["win_id", "chrom", "win_start", "win_end", "win_len"]].merge(
        valid_bp,
        on="win_id",
        how="left",
    )
    validity["valid_bp"] = validity["valid_bp"].fillna(0).astype(int)
  
    # <= длины окна
    validity["valid_bp"] = validity[["valid_bp", "win_len"]].min(axis=1)
    validity["valid_frac"] = validity["valid_bp"] / validity["win_len"]
    validity["is_valid"] = (validity["valid_frac"] >= min_valid_frac).astype(int)
    return validity

def merge_into_full(
    windows_full_tsv: Path,
    validity_df: pd.DataFrame,
    out_path: Path,
) -> pd.DataFrame:
    full_df = pd.read_csv(windows_full_tsv, sep="\t", low_memory=False)
    if "win_id" not in full_df.columns:
        raise RuntimeError(f"'win_id' not found in {windows_full_tsv}")

    merged = full_df.merge(
        validity_df[["win_id", "valid_bp", "valid_frac", "is_valid"]],
        on="win_id",
        how="left",
    )
    merged["valid_bp"] = merged["valid_bp"].fillna(0).astype(int)
    merged["valid_frac"] = merged["valid_frac"].fillna(0.0)
    merged["is_valid"] = merged["is_valid"].fillna(0).astype(int)
    merged.to_csv(out_path, sep="\t", index=False)
    return merged

def main() -> int:
    args = parse_args()
    check_dependencies()
    chrom = str(args.chrom)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    windows_bed = Path(args.windows_bed)
    callable_bed = Path(args.callable_bed)
    windows_full_tsv = Path(args.windows_full_tsv)

    intersect_out = outdir / f"chr{chrom}_windows_x_callable.bed"
    validity_tsv = outdir / f"chr{chrom}_windows_validity.tsv"
    full_valid_tsv = outdir / f"chr{chrom}_windows_full_valid.tsv"

    print(f"[INFO] Chromosome: {chrom}")
    print(f"[INFO] Windows BED: {windows_bed}")
    print(f"[INFO] Callable BED: {callable_bed}")
    print(f"[INFO] Input full TSV: {windows_full_tsv}")
    print(f"[INFO] min_valid_frac: {args.min_valid_frac}")

    windows_df = load_windows(windows_bed, chrom)
    print(f"[INFO] Total windows: {len(windows_df)}")
    run_intersect(windows_bed, callable_bed, intersect_out)
    print(f"[INFO] Saved: {intersect_out}")

    validity_df = compute_validity(
        windows_df=windows_df,
        intersect_path=intersect_out,
        min_valid_frac=args.min_valid_frac,
    )
    validity_df.to_csv(validity_tsv, sep="\t", index=False)
    print(f"[INFO] Saved: {validity_tsv}")
  
    merged = merge_into_full(
        windows_full_tsv=windows_full_tsv,
        validity_df=validity_df,
        out_path=full_valid_tsv,
    )
    print(f"[INFO] Saved: {full_valid_tsv}")
    print("[INFO] Validity summary:")
    print(merged["is_valid"].value_counts().sort_index().to_string())
    print("[INFO] valid_frac summary:")
    print(merged["valid_frac"].describe().to_string())
    return 0

if __name__ == "__main__":
    sys.exit(main())
