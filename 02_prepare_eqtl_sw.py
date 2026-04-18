#!/usr/bin/env python3
# variant_id, gene_id, slope, slope_se, pval_nominal

from __future__ import annotations

import argparse
import glob
import os
import re
import shutil
import subprocess
import sys
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

import numpy as np
import pandas as pd

EMPTY_EQTL_COLUMNS = [
    "variant_id", "gene_id", "slope", "slope_se", "pval_nominal", "tissue",
    "vchr", "vpos", "vref", "valt", "vbuild", "vchr_norm"
]

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Prepare GTEx eQTL and compute Sw per window/tissue."
    )
    parser.add_argument(
        "--chr",
        required=True,
        dest="chrom",
        help="Chromosome name, e.g. 6 or X",
    )
    parser.add_argument(
        "--gtex-dir",
        required=True,
        help="Directory with GTEx *.signif_variant_gene_pairs files",
    )
    parser.add_argument(
        "--windows-bed",
        required=True,
        help="BED with windows: chrom, start, end, win_id",
    )
    parser.add_argument(
        "--fw-tsv",
        required=True,
        help="TSV from step 1 with window-level Fw",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory",
    )
    parser.add_argument(
        "--sw-formula",
        default="max_absz",
        choices=["max_absz", "mean_absz", "sum_absz", "max_z", "mean_z"],
        help="Formula for Sw aggregation (default: max_absz)",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of worker processes for per-tissue GTEx parsing (default: 1)",
    )
    parser.add_argument(
        "--keep-all-eqtl",
        action="store_true",
        help="Keep and save full filtered chromosome-specific eQTL table",
    )
    return parser.parse_args()

def check_dependencies() -> None:
    if shutil.which("bedtools") is None:
        raise RuntimeError("bedtools not found in PATH")

def infer_tissue_name(filepath: str) -> str:
    name = os.path.basename(filepath)
    suffixes = [
        ".v7.signif_variant_gene_pairs.txt.gz",
        ".v7.signif_variant_gene_pairs.txt",
        ".signif_variant_gene_pairs.txt.gz",
        ".signif_variant_gene_pairs.txt",
    ]
    for suffix in suffixes:
        if name.endswith(suffix):
            return name[: -len(suffix)]
    return name

def find_eqtl_files(gtex_dir: Path) -> list[str]:
    patterns = [
        str(gtex_dir / "*.signif_variant_gene_pairs.txt.gz"),
        str(gtex_dir / "*.signif_variant_gene_pairs.txt"),
    ]
    files: list[str] = []
    for pattern in patterns:
        files.extend(glob.glob(pattern))
    files = sorted(set(files))
    if not files:
        raise FileNotFoundError(f"No GTEx signif_variant_gene_pairs files found in {gtex_dir}")
    return files

def normalize_chr_labels(series: pd.Series) -> pd.Series:
    s = series.astype(str)
    return s.str.replace(r"^chr", "", regex=True)

def parse_variant_id(df: pd.DataFrame) -> pd.DataFrame:
    split_df = df["variant_id"].astype(str).str.split("_", n=4, expand=True)
    if split_df.shape[1] < 5:
        raise ValueError("variant_id does not have expected format: chrom_pos_ref_alt_build")
    split_df.columns = ["vchr", "vpos", "vref", "valt", "vbuild"]
    out = pd.concat([df.reset_index(drop=True), split_df.reset_index(drop=True)], axis=1)
    out["vchr"] = out["vchr"].astype(str)
    out["vpos"] = pd.to_numeric(out["vpos"], errors="coerce")
    return out

def process_one_tissue_file(task: tuple[str, str]) -> pd.DataFrame:
    f, chrom = task
    tissue = infer_tissue_name(f)
    compression = "gzip" if f.endswith(".gz") else None
    required_cols = {"variant_id", "gene_id", "slope", "slope_se", "pval_nominal"}
    df = pd.read_csv(f, sep="\t", compression=compression)
    missing = required_cols - set(df.columns)
    if missing:
        raise ValueError(f"{f} is missing required columns: {sorted(missing)}")

    df = df[list(required_cols)].copy()
    df["tissue"] = tissue
    df = parse_variant_id(df)
    df["vchr_norm"] = normalize_chr_labels(df["vchr"])
    df = df[df["vchr_norm"] == str(chrom)].copy()
    if df.empty:
        return pd.DataFrame(columns=EMPTY_EQTL_COLUMNS)

    df = df[
        (df["vref"].str.len() == 1) &
        (df["valt"].str.len() == 1) &
        (df["vref"].str.fullmatch(r"[ACGT]")) &
        (df["valt"].str.fullmatch(r"[ACGT]"))
    ].copy()
    if df.empty:
        return pd.DataFrame(columns=EMPTY_EQTL_COLUMNS)

    df["vpos"] = pd.to_numeric(df["vpos"], errors="coerce")
    df["slope"] = pd.to_numeric(df["slope"], errors="coerce")
    df["slope_se"] = pd.to_numeric(df["slope_se"], errors="coerce")
    df["pval_nominal"] = pd.to_numeric(df["pval_nominal"], errors="coerce")
    df = df.dropna(subset=["vpos", "slope", "slope_se", "pval_nominal"]).copy()
    df = df[df["vpos"] >= 1].copy()

    if df.empty:
        return pd.DataFrame(columns=EMPTY_EQTL_COLUMNS)
    return df

def read_and_filter_gtex(files: list[str], chrom: str, threads: int = 1) -> pd.DataFrame:
    tasks = [(f, chrom) for f in files]
    if threads <= 1:
        kept = [process_one_tissue_file(task) for task in tasks]
    else:
        with ProcessPoolExecutor(max_workers=threads) as ex:
            kept = list(ex.map(process_one_tissue_file, tasks))

    kept = [df for df in kept if not df.empty]
    if not kept:
        return pd.DataFrame(columns=EMPTY_EQTL_COLUMNS)

    eqtl = pd.concat(kept, ignore_index=True)
    return eqtl

def prepare_eqtl_bed(eqtl: pd.DataFrame, chrom: str, eqtl_bed: Path) -> None:
    if eqtl.empty:
        eqtl_bed.write_text("", encoding="utf-8")
        return

    bed = eqtl.copy()
    bed["bed_chr"] = str(chrom)
    bed["bed_start"] = bed["vpos"].astype(int) - 1
    bed["bed_end"] = bed["vpos"].astype(int)
    out_cols = [
        "bed_chr", "bed_start", "bed_end",
        "variant_id", "gene_id", "tissue",
        "slope", "slope_se", "pval_nominal"
    ]
    bed[out_cols].to_csv(eqtl_bed, sep="\t", index=False, header=False)

def run_bedtools_intersect(eqtl_bed: Path, windows_bed: Path, out_path: Path) -> None:
    if eqtl_bed.stat().st_size == 0:
        out_path.write_text("", encoding="utf-8")
        return

    cmd = [
        "bedtools", "intersect",
        "-a", str(eqtl_bed),
        "-b", str(windows_bed),
        "-wa", "-wb",
    ]
    with out_path.open("w", encoding="utf-8") as out:
        subprocess.run(cmd, stdout=out, check=True)

def read_intersection(intersection_path: Path) -> pd.DataFrame:
    cols = [
        "chr_eqtl", "start_eqtl", "end_eqtl",
        "variant_id", "gene_id", "tissue",
        "slope", "slope_se", "pval_nominal",
        "chr_win", "win_start", "win_end", "win_id",
    ]
    if intersection_path.stat().st_size == 0:
        return pd.DataFrame(columns=cols)

    df = pd.read_csv(
        intersection_path,
        sep="\t",
        header=None,
        names=cols,
        dtype={
            "chr_eqtl": str,
            "start_eqtl": int,
            "end_eqtl": int,
            "variant_id": str,
            "gene_id": str,
            "tissue": str,
            "slope": float,
            "slope_se": float,
            "pval_nominal": float,
            "chr_win": str,
            "win_start": int,
            "win_end": int,
            "win_id": int,
        },
    )
    return df

def proxy_clump_min_pvalue(df: pd.DataFrame) -> pd.DataFrame:
    if df.empty:
        return df.copy()

    group_cols = ["win_id", "tissue", "gene_id"]
    tmp = df.copy()
    tmp["abs_slope"] = tmp["slope"].abs()
    tmp = tmp.sort_values(
        by=["win_id", "tissue", "gene_id", "pval_nominal", "abs_slope", "variant_id"],
        ascending=[True, True, True, True, False, True],
        kind="mergesort",
    )
    lead = tmp.groupby(group_cols, as_index=False).first()
    lead = lead.drop(columns=["abs_slope"])
    return lead

def compute_zscore(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["Z"] = np.where(out["slope_se"] != 0, out["slope"] / out["slope_se"], np.nan)
    out["absZ"] = out["Z"].abs()
    out = out.dropna(subset=["Z"]).copy()
    return out

def aggregate_sw(df: pd.DataFrame, formula: str) -> pd.DataFrame:
    if df.empty:
        return pd.DataFrame(columns=["win_id", "tissue", "Sw"])

    if formula == "max_absz":
        sw = (
            df.groupby(["win_id", "tissue"], as_index=False)["absZ"]
            .max()
            .rename(columns={"absZ": "Sw"})
        )
    elif formula == "mean_absz":
        sw = (
            df.groupby(["win_id", "tissue"], as_index=False)["absZ"]
            .mean()
            .rename(columns={"absZ": "Sw"})
        )
    elif formula == "sum_absz":
        sw = (
            df.groupby(["win_id", "tissue"], as_index=False)["absZ"]
            .sum()
            .rename(columns={"absZ": "Sw"})
        )
    elif formula == "max_z":
        sw = (
            df.groupby(["win_id", "tissue"], as_index=False)["Z"]
            .max()
            .rename(columns={"Z": "Sw"})
        )
    elif formula == "mean_z":
        sw = (
            df.groupby(["win_id", "tissue"], as_index=False)["Z"]
            .mean()
            .rename(columns={"Z": "Sw"})
        )
    else:
        raise ValueError(f"Unsupported --sw-formula: {formula}")
    return sw

def compute_window_level_counts(lead_df: pd.DataFrame) -> pd.DataFrame:
    if lead_df.empty:
        return pd.DataFrame(columns=["win_id", "n_eqtl_lead", "n_genes_lead", "n_tissues_lead"])

    counts = lead_df.groupby("win_id", as_index=False).agg(
        n_eqtl_lead=("variant_id", "size"),
        n_genes_lead=("gene_id", pd.Series.nunique),
        n_tissues_lead=("tissue", pd.Series.nunique),
    )
    return counts

def pivot_sw_long_to_wide(sw_long: pd.DataFrame) -> pd.DataFrame:
    if sw_long.empty:
        return pd.DataFrame(columns=["win_id"])

    sw_wide = sw_long.pivot(index="win_id", columns="tissue", values="Sw").reset_index()
    sw_wide.columns.name = None
    return sw_wide

def merge_with_fw(
    fw_tsv: Path,
    sw_long: pd.DataFrame,
    sw_wide: pd.DataFrame,
    count_df: pd.DataFrame,
) -> pd.DataFrame:
    fw = pd.read_csv(fw_tsv, sep="\t")

    if sw_long.empty:
        fw["Sw_global"] = 0.0
    else:
        sw_global = (
            sw_long.groupby("win_id", as_index=False)["Sw"]
            .max()
            .rename(columns={"Sw": "Sw_global"})
        )
        fw = fw.merge(sw_global, on="win_id", how="left")
        fw["Sw_global"] = fw["Sw_global"].fillna(0.0)
      
    if not sw_wide.empty:
        fw = fw.merge(sw_wide, on="win_id", how="left")
        tissue_cols = [c for c in sw_wide.columns if c != "win_id"]
        fw[tissue_cols] = fw[tissue_cols].fillna(0.0)
    if count_df.empty:
        fw["n_eqtl_lead"] = 0
        fw["n_genes_lead"] = 0
        fw["n_tissues_lead"] = 0
    else:
        fw = fw.merge(count_df, on="win_id", how="left")
        for col in ["n_eqtl_lead", "n_genes_lead", "n_tissues_lead"]:
            fw[col] = fw[col].fillna(0).astype(int)
    return fw

def main() -> int:
    args = parse_args()
    check_dependencies()
    chrom = str(args.chrom)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    gtex_dir = Path(args.gtex_dir)
    windows_bed = Path(args.windows_bed)
    fw_tsv = Path(args.fw_tsv)

    eqtl_all_tsv = outdir / f"chr{chrom}_eqtl_all.tsv"
    eqtl_bed = outdir / f"chr{chrom}_eqtl.bed"
    intersection_tsv = outdir / f"chr{chrom}_eqtl_x_windows.tsv"
    lead_eqtl_tsv = outdir / f"chr{chrom}_lead_eqtl.tsv"
    sw_long_tsv = outdir / f"chr{chrom}_Sw_long.tsv"
    sw_wide_tsv = outdir / f"chr{chrom}_Sw_per_tissue.tsv"
    fw_sw_tsv = outdir / f"chr{chrom}_windows_Fw_Sw.tsv"

    print(f"[INFO] Chromosome: {chrom}")
    print(f"[INFO] GTEx dir: {gtex_dir}")
    print(f"[INFO] Windows BED: {windows_bed}")
    print(f"[INFO] Fw TSV: {fw_tsv}")
    print(f"[INFO] Sw formula: {args.sw_formula}")
    print(f"[INFO] Threads: {args.threads}")

    files = find_eqtl_files(gtex_dir)
    print(f"[INFO] Found GTEx files: {len(files)}")
    eqtl = read_and_filter_gtex(files, chrom, threads=args.threads)
    print(f"[INFO] Filtered eQTL rows on chr{chrom}: {len(eqtl)}")

    if args.keep_all_eqtl:
        eqtl.to_csv(eqtl_all_tsv, sep="\t", index=False)
        print(f"[INFO] Saved: {eqtl_all_tsv}")

    prepare_eqtl_bed(eqtl, chrom, eqtl_bed)
    print(f"[INFO] Saved BED: {eqtl_bed}")
    run_bedtools_intersect(eqtl_bed, windows_bed, intersection_tsv)

    intersect_df = read_intersection(intersection_tsv)
    print(f"[INFO] Intersections eQTL x windows: {len(intersect_df)}")
    intersect_df.to_csv(intersection_tsv, sep="\t", index=False)

    lead_df = proxy_clump_min_pvalue(intersect_df)
    print(f"[INFO] Lead SNP rows after clumping: {len(lead_df)}")

    lead_df = compute_zscore(lead_df)
    print(f"[INFO] Lead SNP rows after Z-score filter: {len(lead_df)}")

    lead_df.to_csv(lead_eqtl_tsv, sep="\t", index=False)
    print(f"[INFO] Saved: {lead_eqtl_tsv}")

    sw_long = aggregate_sw(lead_df, args.sw_formula)
    sw_long.to_csv(sw_long_tsv, sep="\t", index=False)
    print(f"[INFO] Saved: {sw_long_tsv}")

    sw_wide = pivot_sw_long_to_wide(sw_long)
    sw_wide.to_csv(sw_wide_tsv, sep="\t", index=False)
    print(f"[INFO] Saved: {sw_wide_tsv}")

    count_df = compute_window_level_counts(lead_df)
    fw_sw = merge_with_fw(fw_tsv, sw_long, sw_wide, count_df)
    fw_sw.to_csv(fw_sw_tsv, sep="\t", index=False)
    print(f"[INFO] Saved: {fw_sw_tsv}")

    if not sw_long.empty:
        print("[INFO] Sw summary:")
        print(sw_long["Sw"].describe().to_string())
        print("[INFO] n_eqtl_lead summary:")
        print(fw_sw["n_eqtl_lead"].describe().to_string())
    else:
        print("[INFO] No Sw values computed; output tables contain zeros where applicable.")
    return 0

if __name__ == "__main__":
    sys.exit(main())
