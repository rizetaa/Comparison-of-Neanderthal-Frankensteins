#!/usr/bin/env python3

from __future__ import annotations
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from scipy.stats import mannwhitneyu

import argparse
import sys
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

CAT_ORDER = ["Promoter", "Near", "Distal"]

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Stratified control with Mann-Whitney and block bootstrap."
    )
    parser.add_argument("--input", required=True, help="Path to chr*_windows_full_valid.tsv or genome_windows_full_valid.tsv")
    parser.add_argument("--outdir", required=True, help="Output directory")
    parser.add_argument(
        "--plot-type",
        default="violin",
        choices=["box", "violin"],
        help="Plot type",
    )
    parser.add_argument(
        "--median-seg-len-file",
        required=True,
        help="Path to median segment length file",
    )
    parser.add_argument(
        "--n-bootstrap",
        type=int,
        default=10000,
        help="Number of bootstrap iterations",
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of worker processes for bootstrap",
    )
    parser.add_argument(
        "--seed",
        type=int,
        default=42,
        help="Random seed",
    )
    return parser.parse_args()

def load_data(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    required = {
        "chrom", "win_start", "win_end", "win_id",
        "Fw", "Sw_global", "DTSS_category", "seg_id", "is_valid"
    }
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(f"Missing required columns in {path}: {sorted(missing)}")
    return df

def prepare_data(df: pd.DataFrame, median_seg_len_bp: int) -> pd.DataFrame:
    df = df[df["is_valid"] == 1].copy()
    df = df[df["DTSS_category"].isin(CAT_ORDER)].copy()
    df["DTSS_category"] = pd.Categorical(df["DTSS_category"], categories=CAT_ORDER, ordered=True)
    df["introgression_group"] = np.where(df["Fw"] > 0, "Fw>0", "Fw=0")
    df["control_block_id"] = (
        df["chrom"].astype(str)
        + ":"
        + (df["win_start"] // int(median_seg_len_bp)).astype(str)
    )
    return df

def observed_rank_diff(x: np.ndarray, y: np.ndarray) -> float:
    combined = np.concatenate([x, y])
    ranks = pd.Series(combined).rank(method="average").to_numpy()
    rx = ranks[: len(x)]
    ry = ranks[len(x):]
    return float(rx.mean() - ry.mean())

def mann_whitney_stats(x: pd.Series, y: pd.Series) -> tuple[float, float]:
    if len(x) == 0 or len(y) == 0:
        return np.nan, np.nan
    u, p = mannwhitneyu(x, y, alternative="two-sided")
    return float(u), float(p)

def resample_blocks(block_map: dict[str, np.ndarray], target_n: int, rng: np.random.Generator) -> np.ndarray:
    block_ids = list(block_map.keys())
    if not block_ids:
        return np.array([], dtype=float)

    sampled = []
    total = 0
    while total < target_n:
        bid = block_ids[rng.integers(0, len(block_ids))]
        arr = block_map[bid]
        sampled.append(arr)
        total += len(arr)
      
    out = np.concatenate(sampled)
    if len(out) > target_n:
        out = out[:target_n]
    return out

def bootstrap_worker(task):
    intro_block_map, ctrl_block_map, n_intro, n_ctrl, n_iter, seed = task
    rng = np.random.default_rng(seed)
    out = np.empty(n_iter, dtype=float)
    for i in range(n_iter):
        bx = resample_blocks(intro_block_map, n_intro, rng)
        by = resample_blocks(ctrl_block_map, n_ctrl, rng)
        out[i] = observed_rank_diff(bx, by)
    return out

def bootstrap_pvalue(
    intro_vals: np.ndarray,
    intro_blocks: pd.Series,
    ctrl_vals: np.ndarray,
    ctrl_blocks: pd.Series,
    n_bootstrap: int,
    seed: int,
    threads: int = 1,
) -> tuple[float, float]:
    obs = observed_rank_diff(intro_vals, ctrl_vals)

    intro_block_map = {
        bid: grp["Sw_global"].to_numpy()
        for bid, grp in pd.DataFrame({"block": intro_blocks, "Sw_global": intro_vals}).groupby("block")
    }
    ctrl_block_map = {
        bid: grp["Sw_global"].to_numpy()
        for bid, grp in pd.DataFrame({"block": ctrl_blocks, "Sw_global": ctrl_vals}).groupby("block")
    }
    n_intro = len(intro_vals)
    n_ctrl = len(ctrl_vals)
  
    if threads <= 1:
        rng = np.random.default_rng(seed)
        boot_stats = np.empty(n_bootstrap, dtype=float)
        for i in range(n_bootstrap):
            bx = resample_blocks(intro_block_map, n_intro, rng)
            by = resample_blocks(ctrl_block_map, n_ctrl, rng)
            boot_stats[i] = observed_rank_diff(bx, by)
    else:
        chunk_sizes = [n_bootstrap // threads] * threads
        for i in range(n_bootstrap % threads):
            chunk_sizes[i] += 1
        tasks = []
        for j, chunk_n in enumerate(chunk_sizes):
            if chunk_n == 0:
                continue
            tasks.append((
                intro_block_map,
                ctrl_block_map,
                n_intro,
                n_ctrl,
                chunk_n,
                seed + j + 1,
            ))
        with ProcessPoolExecutor(max_workers=threads) as ex:
            parts = list(ex.map(bootstrap_worker, tasks))
        boot_stats = np.concatenate(parts)
    p_boot = float((np.sum(np.abs(boot_stats) >= abs(obs)) + 1) / (n_bootstrap + 1))
    return obs, p_boot

def compute_stats(df: pd.DataFrame, n_bootstrap: int, seed: int, threads: int) -> pd.DataFrame:
    rows = []
    for i, cat in enumerate(CAT_ORDER):
        sub = df[df["DTSS_category"] == cat].copy()
        intro = sub[sub["introgression_group"] == "Fw>0"].copy()
        ctrl = sub[sub["introgression_group"] == "Fw=0"].copy()

        x = intro["Sw_global"].dropna()
        y = ctrl["Sw_global"].dropna()
        u_stat, p_mw = mann_whitney_stats(x, y)
      
        if len(x) == 0 or len(y) == 0:
            rank_diff = np.nan
            p_boot = np.nan
        else:
            rank_diff, p_boot = bootstrap_pvalue(
                intro_vals=x.to_numpy(),
                intro_blocks=intro.loc[x.index, "seg_id"],
                ctrl_vals=y.to_numpy(),
                ctrl_blocks=ctrl.loc[y.index, "control_block_id"],
                n_bootstrap=n_bootstrap,
                seed=seed + i,
                threads=threads,
            )
        rows.append({
            "DTSS_category": cat,
            "n_Fw0": int(len(y)),
            "n_Fw_gt_0": int(len(x)),
            "median_Fw0": float(y.median()) if len(y) else np.nan,
            "median_Fw_gt_0": float(x.median()) if len(x) else np.nan,
            "mean_Fw0": float(y.mean()) if len(y) else np.nan,
            "mean_Fw_gt_0": float(x.mean()) if len(x) else np.nan,
            "mannwhitney_u": u_stat,
            "pvalue_mw": p_mw,
            "rank_diff_intro_minus_ctrl": rank_diff,
            "pvalue_bootstrap": p_boot,
        })
    return pd.DataFrame(rows)

def p_to_stars(p: float) -> str:
    if pd.isna(p):
        return "ns"
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "ns"

def save_qc(df: pd.DataFrame, stats_df: pd.DataFrame, out_path: Path, median_seg_len_bp: int, n_bootstrap: int, threads: int) -> None:
    lines = []
    lines.append(f"n_total_windows\t{len(df)}")
    lines.append(f"median_seg_len_bp\t{median_seg_len_bp}")
    lines.append(f"n_bootstrap\t{n_bootstrap}")
    lines.append(f"threads\t{threads}")
    lines.append("")
    lines.append("DTSS_category_counts")
    lines.append(df["DTSS_category"].value_counts().reindex(CAT_ORDER).to_string())
    lines.append("")
    lines.append("introgression_group_counts")
    lines.append(df["introgression_group"].value_counts().to_string())
    lines.append("")
    lines.append("stratified_stats")
    lines.append(stats_df.to_string(index=False))
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

def make_plot(df: pd.DataFrame, stats_df: pd.DataFrame, out_path: Path, plot_type: str) -> None:
    positions = []
    data = []
    labels = []
    base_positions = [1, 4, 7]
    offset = 0.55
    for base, cat in zip(base_positions, CAT_ORDER):
        vals0 = df.loc[
            (df["DTSS_category"] == cat) & (df["introgression_group"] == "Fw=0"),
            "Sw_global"
        ].values
        vals1 = df.loc[
            (df["DTSS_category"] == cat) & (df["introgression_group"] == "Fw>0"),
            "Sw_global"
        ].values
        positions.extend([base - offset, base + offset])
        data.extend([vals0, vals1])
        labels.extend([f"{cat}\nFw=0", f"{cat}\nFw>0"])

    plt.figure(figsize=(10, 5))
    if plot_type == "violin":
        plt.violinplot(data, positions=positions, showmedians=True)
        plt.xticks(positions, labels, rotation=20)
    else:
        plt.boxplot(data, positions=positions, tick_labels=labels)
        plt.xticks(rotation=20)
    for base, cat in zip(base_positions, CAT_ORDER):
        row = stats_df.loc[stats_df["DTSS_category"] == cat].iloc[0]
        star = p_to_stars(row["pvalue_bootstrap"])
        ymax = df.loc[df["DTSS_category"] == cat, "Sw_global"].quantile(0.98)
        plt.text(base, ymax, star, ha="center", va="bottom", fontsize=12)
    plt.ylabel("Sw")
    plt.title("Stratified control: Fw=0 vs Fw>0 within DTSS categories")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()

def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    median_seg_len_bp = int(Path(args.median_seg_len_file).read_text().strip())
    df = load_data(input_path)
    df = prepare_data(df, median_seg_len_bp=median_seg_len_bp)
    stats_df = compute_stats(
        df,
        n_bootstrap=args.n_bootstrap,
        seed=args.seed,
        threads=args.threads,
    )
    stats_df.to_csv(outdir / "stratified_control_stats.tsv", sep="\t", index=False)
    df.to_csv(outdir / "stratified_control_input.tsv", sep="\t", index=False)
    save_qc(
        df,
        stats_df,
        outdir / "stratified_control_qc.txt",
        median_seg_len_bp=median_seg_len_bp,
        n_bootstrap=args.n_bootstrap,
        threads=args.threads,
    )
    make_plot(df, stats_df, outdir / f"stratified_control_{args.plot_type}.png", args.plot_type)
    print(stats_df.to_string(index=False))
    return 0

if __name__ == "__main__":
    sys.exit(main())
