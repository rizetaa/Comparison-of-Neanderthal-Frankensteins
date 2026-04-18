#!/usr/bin/env python3

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Extract TSS from GTF and compute DTSS per window."
    )
    parser.add_argument(
        "--chr",
        required=True,
        dest="chrom",
        help="Chromosome name, e.g. 6 or X",
    )
    parser.add_argument(
        "--gtf",
        required=True,
        help="Path to GENCODE GTF file",
    )
    parser.add_argument(
        "--windows-bed",
        required=True,
        help="Path to windows BED: chrom, start, end, win_id",
    )
    parser.add_argument(
        "--fw-sw-tsv",
        required=True,
        help="Path to table from step 2 (windows_Fw_Sw.tsv)",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory",
    )
    return parser.parse_args()

def check_dependencies() -> None:
    if shutil.which("bedtools") is None:
        raise RuntimeError("bedtools not found in PATH")

def normalize_chr_label(value: str) -> str:
    value = str(value)
    return re.sub(r"^chr", "", value)

def extract_attribute(attr_string: str, key: str) -> str | None:
    pattern = rf'{re.escape(key)} "([^"]+)"'
    match = re.search(pattern, attr_string)
    if match:
        return match.group(1)
    return None

def extract_tss_from_gtf(gtf_path: Path, chrom: str) -> pd.DataFrame:
    chrom = str(chrom)
    rows: list[tuple[str, int, int, str, str | None, str]] = []
    with gtf_path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
              
            seqname, source, feature, start, end, score, strand, frame, attributes = parts
            if feature != "gene":
                continue
            if normalize_chr_label(seqname) != chrom:
                continue
            try:
                start_i = int(start)
                end_i = int(end)
            except ValueError:
                continue
            if strand not in {"+", "-"}:
                continue

            gene_id = extract_attribute(attributes, "gene_id")
            gene_name = extract_attribute(attributes, "gene_name")
            if gene_id is None:
                continue

            # GTF, 1-based
            # BED, 0-based
            # TSS, [p-1, p)
            tss_pos = start_i if strand == "+" else end_i
            bed_start = tss_pos - 1
            bed_end = tss_pos
            if bed_start < 0:
                continue

            rows.append((chrom, bed_start, bed_end, gene_id, gene_name, strand))

    if not rows:
        raise RuntimeError(f"No gene TSS records found for chromosome {chrom} in {gtf_path}")

    tss_df = pd.DataFrame(
        rows,
        columns=["chrom", "tss_start", "tss_end", "gene_id", "gene_name", "strand"],
    )
    tss_df = (
        tss_df
        .sort_values(["chrom", "tss_start", "tss_end", "gene_id"], kind="mergesort")
        .reset_index(drop=True)
    )
    return tss_df

def load_and_sort_windows(windows_bed: Path, chrom: str) -> pd.DataFrame:
    windows = pd.read_csv(
        windows_bed,
        sep="\t",
        header=None,
        names=["chrom", "win_start", "win_end", "win_id"],
        dtype={"chrom": str, "win_start": int, "win_end": int, "win_id": int},
    )
    windows["chrom_norm"] = windows["chrom"].map(normalize_chr_label)
    windows = windows.loc[windows["chrom_norm"] == str(chrom)].copy()
    windows = windows.drop(columns=["chrom_norm"])
    if windows.empty:
        raise RuntimeError(f"No windows found for chromosome {chrom} in {windows_bed}")

    windows = (
        windows
        .sort_values(["chrom", "win_start", "win_end", "win_id"], kind="mergesort")
        .reset_index(drop=True)
    )
    return windows

def run_bedtools_closest(windows_sorted_bed: Path, tss_bed: Path, output_path: Path) -> None:
    cmd = [
        "bedtools",
        "closest",
        "-a",
        str(windows_sorted_bed),
        "-b",
        str(tss_bed),
        "-d",
        "-t",
        "first",
    ]
    print(f"[RUN] {' '.join(cmd)}")
    with output_path.open("w", encoding="utf-8") as out:
        subprocess.run(cmd, stdout=out, check=True)

def assign_dtss_category(dtss: int) -> str:
    if dtss < 5000:
        return "Promoter"
    if dtss < 50000:
        return "Near"
    return "Distal"

def main() -> int:
    args = parse_args()
    check_dependencies()
    chrom = str(args.chrom)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    gtf_path = Path(args.gtf)
    windows_bed = Path(args.windows_bed)
    fw_sw_tsv = Path(args.fw_sw_tsv)

    tss_bed = outdir / f"chr{chrom}_tss.bed"
    windows_sorted_bed = outdir / f"chr{chrom}_windows.sorted.bed"
    closest_out = outdir / f"chr{chrom}_windows_x_tss.bed"
    dtss_tsv = outdir / f"chr{chrom}_windows_DTSS.tsv"
    merged_tsv = outdir / f"chr{chrom}_windows_Fw_Sw_DTSS.tsv"

    print(f"[INFO] Chromosome: {chrom}")
    print(f"[INFO] GTF: {gtf_path}")
    print(f"[INFO] Windows BED: {windows_bed}")
    print(f"[INFO] Fw+Sw table: {fw_sw_tsv}")

    # Извлекаем и сохраняем отсорченные TSS
    print("[INFO] Extracting gene TSS from GTF")
    tss_df = extract_tss_from_gtf(gtf_path, chrom)
    tss_df.to_csv(
        tss_bed,
        sep="\t",
        index=False,
        header=False,
        columns=["chrom", "tss_start", "tss_end", "gene_id", "gene_name", "strand"],
    )
    print(f"[INFO] TSS count: {len(tss_df)}")
    print(f"[INFO] Saved: {tss_bed}")

    # Загружаем и сохраняем отсорченные окна
    windows_df = load_and_sort_windows(windows_bed, chrom)
    windows_df.to_csv(
        windows_sorted_bed,
        sep="\t",
        index=False,
        header=False,
        columns=["chrom", "win_start", "win_end", "win_id"],
    )
    print(f"[INFO] Window count: {len(windows_df)}")
    print(f"[INFO] Saved: {windows_sorted_bed}")
    # Ближайший TSS
    run_bedtools_closest(windows_sorted_bed, tss_bed, closest_out)
    print(f"[INFO] Saved: {closest_out}")
  
    closest_cols = [
        "chr_win",
        "win_start",
        "win_end",
        "win_id",
        "chr_tss",
        "tss_start",
        "tss_end",
        "gene_id",
        "gene_name",
        "strand",
        "distance",
    ]
    closest_df = pd.read_csv(
        closest_out,
        sep="\t",
        header=None,
        names=closest_cols,
        dtype={
            "chr_win": str,
            "win_start": int,
            "win_end": int,
            "win_id": int,
            "chr_tss": str,
            "tss_start": int,
            "tss_end": int,
            "gene_id": str,
            "gene_name": str,
            "strand": str,
            "distance": int,
        },
    )
    if len(closest_df) != len(windows_df):
        raise RuntimeError(
            f"Expected one closest-TSS record per window, got {len(closest_df)} "
            f"for {len(windows_df)} windows"
        )
    closest_df["DTSS"] = closest_df["distance"]
    closest_df["DTSS_category"] = closest_df["DTSS"].map(assign_dtss_category)
    dtss_df = closest_df[
        ["win_id", "gene_id", "gene_name", "strand", "DTSS", "DTSS_category"]
    ].copy()
    dtss_df.to_csv(dtss_tsv, sep="\t", index=False)
    print(f"[INFO] Saved: {dtss_tsv}")

    # Объединяем с Fw+Sw
    fw_sw_df = pd.read_csv(fw_sw_tsv, sep="\t")
    if "win_id" not in fw_sw_df.columns:
        raise RuntimeError(f"'win_id' column not found in {fw_sw_tsv}")

    merged_df = fw_sw_df.merge(dtss_df, on="win_id", how="left")
    if merged_df["DTSS"].isna().any():
        missing = int(merged_df["DTSS"].isna().sum())
        raise RuntimeError(f"DTSS missing for {missing} windows after merge")

    merged_df.to_csv(merged_tsv, sep="\t", index=False)
    print(f"[INFO] Saved: {merged_tsv}")
    print("[INFO] DTSS category counts:")
    print(merged_df["DTSS_category"].value_counts().to_string())
    print("[INFO] DTSS summary:")
    print(merged_df["DTSS"].describe().to_string())
    return 0

if __name__ == "__main__":
    sys.exit(main())
