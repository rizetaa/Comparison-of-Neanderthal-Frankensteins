#!/usr/bin/env python3

from __future__ import annotations
from pathlib import Path

import argparse
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

FREQ_ORDER = ["Zero", "Rare", "Low", "Intermediate", "High", "Very_High"]

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="QC for introgression landscape: iSFS-style bin counts and Manhattan plot."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to chr*_windows_full_valid.tsv or genome_windows_full_valid.tsv",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory",
    )
    parser.add_argument(
        "--valid-only",
        action="store_true",
        help="Use only valid windows (is_valid == 1)",
    )
    parser.add_argument(
        "--point-size",
        type=float,
        default=4.0,
        help="Point size for Manhattan plot",
    )
    return parser.parse_args()

def load_data(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    required = {"chrom", "win_start", "win_end", "Fw", "freq_bin"}
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(f"Missing required columns in {path}: {sorted(missing)}")

    if "is_valid" not in df.columns:
        df["is_valid"] = 1
    df["chrom"] = df["chrom"].astype(str)
    return df

def filter_data(df: pd.DataFrame, valid_only: bool) -> pd.DataFrame:
    out = df.copy()
    if valid_only:
        out = out[out["is_valid"] == 1].copy()
    if out.empty:
        raise RuntimeError("No windows left after filtering")
    return out

def save_bin_counts(df: pd.DataFrame, out_path: Path) -> pd.DataFrame:
    counts = (
        df["freq_bin"]
        .value_counts()
        .reindex(FREQ_ORDER, fill_value=0)
        .reset_index()
    )
    counts.columns = ["freq_bin", "n_windows"]
    counts.to_csv(out_path, sep="\t", index=False)
    return counts

def plot_bin_counts(counts: pd.DataFrame, out_path: Path) -> None:
    plt.figure(figsize=(7, 4.5))
    x = np.arange(len(counts))
    plt.bar(x, counts["n_windows"].values)
    plt.xticks(x, counts["freq_bin"].values, rotation=20)
    plt.ylabel("Number of windows")
    plt.xlabel("Introgression frequency bin")
    plt.title("Spectrum of introgression frequency")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()

def make_cumulative_positions(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    chrom_order = sorted(df["chrom"].unique(), key=lambda x: (not x.isdigit(), int(x) if x.isdigit() else x))
    chrom_sizes = (
        df.groupby("chrom", as_index=False)["win_end"]
        .max()
        .rename(columns={"win_end": "chrom_end"})
    )
    chrom_sizes["chrom"] = pd.Categorical(chrom_sizes["chrom"], categories=chrom_order, ordered=True)
    chrom_sizes = chrom_sizes.sort_values("chrom").reset_index(drop=True)
    offsets = []
    cum = 0
    for _, row in chrom_sizes.iterrows():
        offsets.append((str(row["chrom"]), cum))
        cum += int(row["chrom_end"])

    offset_df = pd.DataFrame(offsets, columns=["chrom", "offset"])
    df2 = df.merge(offset_df, on="chrom", how="left")
    df2["plot_x"] = df2["win_start"] + df2["offset"]
    tick_df = (
        df2.groupby("chrom", as_index=False)
        .agg(xmin=("plot_x", "min"), xmax=("plot_x", "max"))
    )
    tick_df["xtick"] = (tick_df["xmin"] + tick_df["xmax"]) / 2
    return df2, tick_df

def plot_manhattan(df: pd.DataFrame, out_path: Path, point_size: float) -> None:
    df2, tick_df = make_cumulative_positions(df)
    plt.figure(figsize=(11, 4.5))
    chrom_order = tick_df["chrom"].tolist()
    for i, chrom in enumerate(chrom_order):
        sub = df2[df2["chrom"] == chrom]
        plt.scatter(
            sub["plot_x"],
            sub["Fw"],
            s=point_size,
            alpha=0.7,
            rasterized=True,
        )
    plt.xticks(tick_df["xtick"], tick_df["chrom"], rotation=0)
    plt.xlabel("Genomic position")
    plt.ylabel("Fw")
    plt.title("Manhattan plot of introgression frequency")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()

def save_qc(df: pd.DataFrame, counts: pd.DataFrame, out_path: Path, valid_only: bool) -> None:
    lines = []
    lines.append(f"n_windows\t{len(df)}")
    lines.append(f"valid_only\t{int(valid_only)}")
    lines.append(f"n_chromosomes\t{df['chrom'].nunique()}")
    lines.append("")
    lines.append("freq_bin_counts")
    lines.append(counts.to_string(index=False))
    lines.append("")
    lines.append("Fw_summary")
    lines.append(df["Fw"].describe().to_string())
    lines.append("")
    lines.append("chrom_counts")
    lines.append(df["chrom"].value_counts().sort_index().to_string())
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    df = load_data(input_path)
    df = filter_data(df, valid_only=args.valid_only)

    counts = save_bin_counts(df, outdir / "introgression_bin_counts.tsv")
    plot_bin_counts(counts, outdir / "introgression_bin_barplot.png")
    plot_manhattan(df, outdir / "introgression_manhattan.png", point_size=args.point_size)
    save_qc(df, counts, outdir / "introgression_qc.txt", valid_only=args.valid_only)
  
    print(f"Windows used: {len(df)}")
    print(f"Chromosomes: {df['chrom'].nunique()}")
    print(counts.to_string(index=False))
    return 0

if __name__ == "__main__":
    sys.exit(main())
