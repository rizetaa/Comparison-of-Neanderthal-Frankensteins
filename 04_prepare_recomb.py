#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Compute mean recombination rate per window from HapMap map."
    )
    parser.add_argument(
        "--chr",
        required=True,
        dest="chrom",
        help="Chromosome name, e.g. 6 or X",
    )
    parser.add_argument(
        "--hapmap",
        required=True,
        help="Path to HapMap recombination map for one chromosome",
    )
    parser.add_argument(
        "--fw-sw-dtss-tsv",
        required=True,
        help="Path to chr*_windows_Fw_Sw_DTSS.tsv",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory",
    )
    return parser.parse_args()

def normalize_chr_label(value: str) -> str:
    return re.sub(r"^chr", "", str(value))

def infer_position_column(columns: list[str]) -> str:
    priority = [
        "position",
        "pos",
        "bp",
        "Position(bp)",
        "position(bp)",
        "Position",
    ]
    colmap = {c.lower(): c for c in columns}
    for candidate in priority:
        if candidate.lower() in colmap:
            return colmap[candidate.lower()]

    for c in columns:
        cl = c.lower()
        if "pos" in cl or "bp" in cl:
            return c

    raise ValueError(
        f"Could not infer position column from HapMap file columns: {columns}"
    )

def infer_rate_column(columns: list[str]) -> str:
    priority = [
        "rate",
        "rate(cM/Mb)",
        "rate(cm/mb)",
        "recombrate",
        "recomb_rate",
    ]
    colmap = {c.lower(): c for c in columns}
    for candidate in priority:
        if candidate.lower() in colmap:
            return colmap[candidate.lower()]

    for c in columns:
        cl = c.lower()
        if "rate" in cl and ("cm" in cl or "mb" in cl or "recomb" in cl):
            return c

    for c in columns:
        if "rate" in c.lower():
            return c

    raise ValueError(
        f"Could not infer recombination-rate column from HapMap file columns: {columns}"
    )

def infer_chrom_column(columns: list[str]) -> str | None:
    priority = ["chr", "chrom", "chromosome"]
    colmap = {c.lower(): c for c in columns}
    for candidate in priority:
        if candidate.lower() in colmap:
            return colmap[candidate.lower()]
    return None

def read_hapmap_map(hapmap_path: Path, chrom: str) -> pd.DataFrame:
    df = pd.read_csv(hapmap_path, sep=r"\s+|,", engine="python")
    if df.empty:
        raise RuntimeError(f"HapMap file is empty: {hapmap_path}")

    pos_col = infer_position_column(list(df.columns))
    rate_col = infer_rate_column(list(df.columns))
    chr_col = infer_chrom_column(list(df.columns))
    keep_cols = [pos_col, rate_col] + ([chr_col] if chr_col else [])
    df = df[keep_cols].copy()

    rename_map = {
        pos_col: "position_bp",
        rate_col: "recomb_rate_cM_Mb",
    }
    if chr_col:
        rename_map[chr_col] = "chrom"

    df = df.rename(columns=rename_map)
    df["position_bp"] = pd.to_numeric(df["position_bp"], errors="coerce")
    df["recomb_rate_cM_Mb"] = pd.to_numeric(df["recomb_rate_cM_Mb"], errors="coerce")
    df = df.dropna(subset=["position_bp", "recomb_rate_cM_Mb"]).copy()
    df["position_bp"] = df["position_bp"].astype(int)

    if chr_col:
        df["chrom"] = df["chrom"].map(normalize_chr_label)
        df = df.loc[df["chrom"] == str(chrom)].copy()
    if df.empty:
        raise RuntimeError(f"No valid HapMap rows found for chromosome {chrom}")

    df = (
        df.sort_values("position_bp", kind="mergesort")
        .drop_duplicates(subset=["position_bp"], keep="first")
        .reset_index(drop=True)
    )
    return df[["position_bp", "recomb_rate_cM_Mb"]]

def load_windows_table(fw_sw_dtss_tsv: Path, chrom: str) -> pd.DataFrame:
    df = pd.read_csv(fw_sw_dtss_tsv, sep="\t")
    required = {"chrom", "win_start", "win_end", "win_id"}
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(
            f"Missing required columns in {fw_sw_dtss_tsv}: {sorted(missing)}"
        )
      
    df["chrom_norm"] = df["chrom"].map(normalize_chr_label)
    df = df.loc[df["chrom_norm"] == str(chrom)].copy()
    df = df.drop(columns=["chrom_norm"])
    if df.empty:
        raise RuntimeError(f"No windows found for chromosome {chrom} in {fw_sw_dtss_tsv}")

    df = (
        df.sort_values(["win_start", "win_end", "win_id"], kind="mergesort")
        .reset_index(drop=True)
    )
    return df

def compute_window_recomb(
    windows_df: pd.DataFrame,
    recomb_df: pd.DataFrame,
) -> pd.DataFrame:
    positions = recomb_df["position_bp"].to_numpy()
    rates = recomb_df["recomb_rate_cM_Mb"].to_numpy()
    mean_rates: list[float] = []
    n_points: list[int] = []

    for start, end in zip(windows_df["win_start"], windows_df["win_end"]):
        # [start, end)
        left = np.searchsorted(positions, start, side="left")
        right = np.searchsorted(positions, end, side="left")
        if right > left:
            vals = rates[left:right]
            mean_rates.append(float(np.mean(vals)))
            n_points.append(int(right - left))
        else:
            mean_rates.append(np.nan)
            n_points.append(0)

    out = windows_df[["win_id", "chrom", "win_start", "win_end"]].copy()
    out["recomb_rate_cM_Mb_raw"] = mean_rates
    out["recomb_n_points"] = n_points
    # Заполняем окна без внутренних точек из ближайших соседних окон
    out["recomb_rate_cM_Mb"] = (
        out["recomb_rate_cM_Mb_raw"]
        .ffill()
        .bfill()
    )
    if out["recomb_rate_cM_Mb"].isna().all():
        raise RuntimeError("All recombination rates are NaN after fill")
    return out

def main() -> int:
    args = parse_args()
    chrom = str(args.chrom)
    hapmap_path = Path(args.hapmap)
    fw_sw_dtss_tsv = Path(args.fw_sw_dtss_tsv)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    recomb_tsv = outdir / f"chr{chrom}_windows_recomb.tsv"
    full_tsv = outdir / f"chr{chrom}_windows_full.tsv"
    print(f"[INFO] Chromosome: {chrom}")
    print(f"[INFO] HapMap file: {hapmap_path}")
    print(f"[INFO] Input table: {fw_sw_dtss_tsv}")

    recomb_df = read_hapmap_map(hapmap_path, chrom)
    print(f"[INFO] HapMap points loaded: {len(recomb_df)}")
    print(
        f"[INFO] HapMap position range: "
        f"{recomb_df['position_bp'].min()} - {recomb_df['position_bp'].max()}"
    )
    windows_df = load_windows_table(fw_sw_dtss_tsv, chrom)
    print(f"[INFO] Windows loaded: {len(windows_df)}")

    recomb_win_df = compute_window_recomb(windows_df, recomb_df)
    recomb_win_df.to_csv(recomb_tsv, sep="\t", index=False)
    print(f"[INFO] Saved: {recomb_tsv}")

    merged_df = windows_df.merge(
        recomb_win_df[["win_id", "recomb_rate_cM_Mb", "recomb_rate_cM_Mb_raw", "recomb_n_points"]],
        on="win_id",
        how="left",
    )
    if merged_df["recomb_rate_cM_Mb"].isna().any():
        missing = int(merged_df["recomb_rate_cM_Mb"].isna().sum())
        raise RuntimeError(f"Missing recombination rate for {missing} windows after merge")

    merged_df.to_csv(full_tsv, sep="\t", index=False)
    print(f"[INFO] Saved: {full_tsv}")
    n_zero_points = int((merged_df["recomb_n_points"] == 0).sum())
    print(f"[INFO] Windows with no internal HapMap points: {n_zero_points}")
    print("[INFO] recomb_rate_cM_Mb summary:")
    print(merged_df["recomb_rate_cM_Mb"].describe().to_string())
    return 0

if __name__ == "__main__":
    sys.exit(main())
