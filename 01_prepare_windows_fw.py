#!/usr/bin/env python3
# Sample, CHROM, Start, End, Length
# Start/End, 1-based

from __future__ import annotations
from pathlib import Path

import argparse
import shutil
import subprocess
import sys
import pandas as pd

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create windows and compute global introgression frequency Fw."
    )
    parser.add_argument(
        "--chr-len-file",
        required=True,
        help="Path to chromosome length file, e.g. hg19.chrom.len",
    )
    parser.add_argument(
        "--chr",
        required=True,
        dest="chrom",
        help="Chromosome name as in hg19.chrom.len and NIS data, e.g. 6 or X",
    )
    parser.add_argument(
        "--nis",
        required=True,
        help="Path to NIS TSV for one chromosome",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory",
    )
    parser.add_argument(
        "--win-size",
        type=int,
        default=1000,
        help="Window size in bp (default: 1000)",
    )
    parser.add_argument(
        "--min-overlap-frac",
        type=float,
        default=0.5,
        help="Minimum fraction of a window covered by a segment (default: 0.5)",
    )
    return parser.parse_args()

def check_dependencies() -> None:
    if shutil.which("bedtools") is None:
        raise RuntimeError("bedtools not found in PATH")

def read_chr_length(chr_len_file: Path, chrom: str) -> int:
    chrom_df = pd.read_csv(
        chr_len_file,
        sep=r"\s+",
        header=None,
        names=["chrom", "length"],
        dtype={"chrom": str, "length": int},
    )
    match = chrom_df.loc[chrom_df["chrom"].astype(str) == str(chrom), "length"]
    if match.empty:
        raise ValueError(f"Chromosome {chrom!r} not found in {chr_len_file}")
    return int(match.iloc[0])

def write_genome_file(genome_file: Path, chrom: str, chrom_len: int) -> None:
    genome_file.write_text(f"{chrom}\t{chrom_len}\n", encoding="utf-8")

def make_windows(genome_file: Path, windows_bed: Path, win_size: int) -> None:
    cmd = [
        "bedtools",
        "makewindows",
        "-g",
        str(genome_file),
        "-w",
        str(win_size),
    ]
    print(f"[RUN] makewindows -> {windows_bed}")
    with windows_bed.open("w", encoding="utf-8") as out:
        subprocess.run(cmd, stdout=out, check=True)

def load_nis(nis_path: Path) -> pd.DataFrame:
    nis = pd.read_csv(nis_path, sep=r"\s+")
    required = {"Sample", "Start", "End"}
    missing = required - set(nis.columns)
    if missing:
        raise ValueError(
            f"NIS file {nis_path} is missing required columns: {sorted(missing)}"
        )
    return nis

def prepare_nis_bed(nis: pd.DataFrame, chrom: str, nis_bed: Path) -> int:
    nis = nis.copy()
    nis["Sample"] = nis["Sample"].astype(str)
    haplotypes = nis["Sample"].dropna().unique()
    H = len(haplotypes)
    if H == 0:
        raise ValueError("No haplotypes found in NIS file")

    bed = pd.DataFrame(
        {
            "chrom": str(chrom),
            "start": nis["Start"].astype(int) - 1,
            "end": nis["End"].astype(int),
            "hap_id": nis["Sample"],
        }
    )
    n_bad = int((bed["start"] < 0).sum())
    if n_bad > 0:
        print(f"[WARN] Found {n_bad} NIS segments with Start < 1; clipping to chromosome start")
        bed.loc[bed["start"] < 0, "start"] = 0
    bed.to_csv(nis_bed, sep="\t", index=False, header=False)
    return H

def save_median_seg_len(nis: pd.DataFrame, outdir: Path, chrom: str) -> int:
    if "Length" not in nis.columns:
        raise ValueError("NIS file does not contain required column 'Length'")
    median_seg_len = int(pd.to_numeric(nis["Length"], errors="coerce").dropna().median())
    (outdir / f"chr{chrom}_median_seg_len.txt").write_text(
        f"{median_seg_len}\n",
        encoding="utf-8"
    )
    return median_seg_len

def intersect_windows_with_nis(
    windows_bed: Path,
    nis_bed: Path,
    intersect_out: Path,
) -> None:
    cmd = [
        "bedtools",
        "intersect",
        "-a",
        str(windows_bed),
        "-b",
        str(nis_bed),
        "-wo",
    ]
    print(f"[RUN] intersect windows x NIS -> {intersect_out}")
    with intersect_out.open("w", encoding="utf-8") as out:
        subprocess.run(cmd, stdout=out, check=True)

def assign_freq_bin(fw: float) -> str:
    if fw == 0:
        return "Zero"
    if fw <= 0.02:
        return "Rare"
    if fw <= 0.05:
        return "Low"
    if fw <= 0.10:
        return "Intermediate"
    if fw <= 0.20:
        return "High"
    return "Very_High"

def compute_fw(
    windows_bed: Path,
    intersect_out: Path,
    H: int,
    win_size: int,
    min_overlap_frac: float,
) -> tuple[pd.DataFrame, int]:
    min_overlap_bp = int(round(win_size * min_overlap_frac))
    windows = pd.read_csv(
        windows_bed,
        sep="\t",
        header=None,
        names=["chrom", "win_start", "win_end"],
        dtype={"chrom": str, "win_start": int, "win_end": int},
    )
    windows["win_id"] = range(len(windows))

    def empty_result() -> pd.DataFrame:
        fw_df = windows.copy()
        fw_df["n_hap_introgressed"] = 0
        fw_df["Fw"] = 0.0
        fw_df["freq_bin"] = "Zero"
        fw_df["seg_id"] = pd.NA
        fw_df["seg_hap_id"] = pd.NA
        fw_df["seg_start"] = pd.NA
        fw_df["seg_end"] = pd.NA
        fw_df["seg_overlap_bp"] = pd.NA
        return fw_df

    if intersect_out.stat().st_size == 0:
        return empty_result(), min_overlap_bp

    cols = [
        "chrom_win",
        "win_start",
        "win_end",
        "chrom_nis",
        "seg_start",
        "seg_end",
        "hap_id",
        "overlap_bp",
    ]
    ix = pd.read_csv(
        intersect_out,
        sep="\t",
        header=None,
        names=cols,
        dtype={
            "chrom_win": str,
            "win_start": int,
            "win_end": int,
            "chrom_nis": str,
            "seg_start": int,
            "seg_end": int,
            "hap_id": str,
            "overlap_bp": int,
        },
    )
    ix = ix.loc[ix["overlap_bp"] >= min_overlap_bp].copy()

    if ix.empty:
        return empty_result(), min_overlap_bp

    ix = ix.merge(
        windows[["chrom", "win_start", "win_end", "win_id"]].rename(
            columns={"chrom": "chrom_win"}
        ),
        on=["chrom_win", "win_start", "win_end"],
        how="left",
        validate="many_to_one",
    )

    if ix["win_id"].isna().any():
        raise RuntimeError("Some intersections could not be matched back to window IDs")

    n_introgressed = (
        ix.groupby("win_id", as_index=False)["hap_id"]
        .nunique()
        .rename(columns={"hap_id": "n_hap_introgressed"})
    )
    ix["seg_id"] = (
        ix["hap_id"].astype(str)
        + ":"
        + ix["seg_start"].astype(str)
        + "-"
        + ix["seg_end"].astype(str)
    )
    rep_seg = (
        ix.sort_values(
            ["win_id", "overlap_bp", "hap_id", "seg_start", "seg_end"],
            ascending=[True, False, True, True, True],
            kind="mergesort",
        )
        .groupby("win_id", as_index=False)
        .first()[["win_id", "seg_id", "hap_id", "seg_start", "seg_end", "overlap_bp"]]
        .rename(
            columns={
                "hap_id": "seg_hap_id",
                "overlap_bp": "seg_overlap_bp",
            }
        )
    )

    fw_df = windows.merge(n_introgressed, on="win_id", how="left")
    fw_df = fw_df.merge(rep_seg, on="win_id", how="left")
    fw_df["n_hap_introgressed"] = fw_df["n_hap_introgressed"].fillna(0).astype(int)
    fw_df["Fw"] = fw_df["n_hap_introgressed"] / H
    fw_df["freq_bin"] = fw_df["Fw"].map(assign_freq_bin)
    return fw_df, min_overlap_bp

def save_outputs(
    fw_df: pd.DataFrame,
    windows_simple_bed: Path,
    fw_tsv: Path,
    h_total_file: Path,
    H: int,
) -> None:
    fw_df.to_csv(fw_tsv, sep="\t", index=False)

    fw_df[["chrom", "win_start", "win_end", "win_id"]].to_csv(
        windows_simple_bed,
        sep="\t",
        index=False,
        header=False,
    )
    h_total_file.write_text(f"{H}\n", encoding="utf-8")


def main() -> int:
    args = parse_args()
    check_dependencies()

    chrom = str(args.chrom)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    chr_len_file = Path(args.chr_len_file)
    nis_path = Path(args.nis)

    genome_file = outdir / f"chr{chrom}.genome"
    windows_bed = outdir / f"chr{chrom}_windows_1kb.bed"
    windows_simple_bed = outdir / f"chr{chrom}_windows.bed"
    nis_bed = outdir / f"chr{chrom}_nis_segments.bed"
    intersect_out = outdir / f"chr{chrom}_windows_x_nis.bed"
    fw_tsv = outdir / f"chr{chrom}_windows_Fw.tsv"
    h_total_file = outdir / f"chr{chrom}_H_total.txt"

    print(f"[INFO] Chromosome: {chrom}")
    print(f"[INFO] Window size: {args.win_size}")
    print(f"[INFO] Output dir: {outdir}")
    chrom_len = read_chr_length(chr_len_file, chrom)
    print(f"[INFO] Chromosome length: {chrom_len}")

    write_genome_file(genome_file, chrom, chrom_len)
    make_windows(genome_file, windows_bed, args.win_size)

    nis = load_nis(nis_path)
    H = prepare_nis_bed(nis, chrom, nis_bed)
    print(f"[INFO] Total haplotypes H = {H}")
    median_seg_len = save_median_seg_len(nis, outdir, chrom)
    #median_seg_len = save_median_seg_len(nis, outdir)
    print(f"[INFO] Median NIS segment length: {median_seg_len}")
    intersect_windows_with_nis(windows_bed, nis_bed, intersect_out)

    fw_df, min_overlap_bp = compute_fw(
        windows_bed=windows_bed,
        intersect_out=intersect_out,
        H=H,
        win_size=args.win_size,
        min_overlap_frac=args.min_overlap_frac,
    )
    save_outputs(
        fw_df=fw_df,
        windows_simple_bed=windows_simple_bed,
        fw_tsv=fw_tsv,
        h_total_file=h_total_file,
        H=H,
    )
    print(f"[INFO] Minimum overlap threshold: {min_overlap_bp} bp")
    print(f"[INFO] Total windows: {len(fw_df)}")
    print("[INFO] Frequency bin counts:")
    print(fw_df["freq_bin"].value_counts().sort_index().to_string())
    print(f"[INFO] Saved: {fw_tsv}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
