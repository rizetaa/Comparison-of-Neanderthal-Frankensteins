#!/usr/bin/env python3

from __future__ import annotations
from pathlib import Path

import argparse
import sys
import matplotlib.pyplot as plt
import pandas as pd

WINDOW_KEY = ["chrom", "win_start", "win_end", "win_id"]
REQUIRED_WINDOWS_COLUMNS = [
    "chrom",
    "win_start",
    "win_end",
    "win_id",
    "Fw",
    "Sw_global",
    "is_valid",
    "n_hap_introgressed",
]
REQUIRED_LEAD_EQTL_COLUMNS = [
    "chrom",
    "win_start",
    "win_end",
    "win_id",
    "gene_id",
    "tissue",
]

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Find adaptive introgression candidate windows and regions.")
    p.add_argument("--input", required=True, help="Genome-wide windows table.")
    p.add_argument("--lead-eqtl", required=True, help="Genome-wide lead-eQTL table.")
    p.add_argument("--outdir", required=True, help="Output directory.")
    p.add_argument("--fw-quantile", type=float, default=0.95, help="Quantile threshold for Fw.")
    p.add_argument("--sw-quantile", type=float, default=0.95, help="Quantile threshold for Sw_global.")
    p.add_argument("--max-gap-bp", type=int, default=0, help="Maximum allowed gap between neighboring candidate windows.")
    return p.parse_args()

def require_columns(df: pd.DataFrame, required: list[str], label: str) -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{label}: missing required columns: {missing}")

def normalize_numeric(df: pd.DataFrame, cols: list[str]) -> pd.DataFrame:
    out = df.copy()
    for col in cols:
        if col in out.columns:
            out[col] = pd.to_numeric(out[col], errors="coerce")
    return out

def chrom_sort_key(series: pd.Series) -> pd.Series:
    x = pd.to_numeric(series, errors="coerce")
    if x.notna().all():
        return x
    return series.astype(str)

def read_windows(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", low_memory=False)
    require_columns(df, REQUIRED_WINDOWS_COLUMNS, "windows")
    numeric_cols = [
        "win_start",
        "win_end",
        "win_id",
        "Fw",
        "Sw_global",
        "is_valid",
        "n_hap_introgressed",
        "seg_start",
        "seg_end",
        "seg_overlap_bp",
        "DTSS",
        "recomb_rate_cM_Mb",
        "recomb_rate_cM_Mb_raw",
        "recomb_n_points",
        "valid_bp",
        "valid_frac",
        "n_eqtl_lead",
        "n_genes_lead",
        "n_tissues_lead",
    ]
    df = normalize_numeric(df, numeric_cols)
    return df

def read_lead_eqtl(path: Path) -> pd.DataFrame:
    # Только нужные столбцы
    all_cols = pd.read_csv(path, sep="\t", nrows=0).columns.tolist()
    wanted = [
        "chrom", "win_start", "win_end", "win_id",
        "gene_id", "tissue",
        "variant_id", "Z", "absZ", "slope", "slope_se",
        "pval_nominal", "chr_eqtl", "start_eqtl", "end_eqtl",
    ]
    usecols = [c for c in wanted if c in all_cols]
    df = pd.read_csv(path, sep="\t", usecols=usecols, low_memory=False)
    require_columns(df, REQUIRED_LEAD_EQTL_COLUMNS, "lead-eQTL")

    numeric_cols = [
        "win_start",
        "win_end",
        "win_id",
        "Z",
        "absZ",
        "slope",
        "slope_se",
        "pval_nominal",
        "start_eqtl",
        "end_eqtl",
    ]
    df = normalize_numeric(df, numeric_cols)
    return df

def uniq_join(series: pd.Series) -> str:
    vals = []
    seen = set()
    for x in series.dropna():
        x = str(x).strip()
        if not x:
            continue
        if x not in seen:
            seen.add(x)
            vals.append(x)
    return ",".join(vals)

def union_gene_lists(series: pd.Series) -> str:
    vals = []
    seen = set()

    for x in series.dropna():
        for g in str(x).split(","):
            g = g.strip()
            if not g:
                continue
            if g not in seen:
                seen.add(g)
                vals.append(g)
    return ",".join(vals)

def summarize_candidate_lead_eqtl(lead_eqtl_subset: pd.DataFrame) -> pd.DataFrame:
    if lead_eqtl_subset.empty:
        return pd.DataFrame(columns=WINDOW_KEY + [
            "n_genes_lead", "target_genes", "n_target_genes",
            "n_tissues_lead", "n_eqtl_lead"
        ])

    x = lead_eqtl_subset.copy()
    if "absZ" in x.columns:
        top_idx = x.groupby(WINDOW_KEY)["absZ"].idxmax()
        top = x.loc[top_idx].copy()
    else:
        top = x.drop_duplicates(WINDOW_KEY).copy()

    keep_cols = []
    for col in [
        "gene_id",
        "tissue",
        "variant_id",
        "Z",
        "absZ",
        "slope",
        "slope_se",
        "pval_nominal",
        "chr_eqtl",
        "start_eqtl",
        "end_eqtl",
    ]:
        if col in top.columns:
            keep_cols.append(col)

    top = top[WINDOW_KEY + keep_cols].copy()
    top = top.rename(columns={c: f"top_{c}" for c in keep_cols})
    genes = (
        x.groupby(WINDOW_KEY)["gene_id"]
        .agg(
            n_genes_lead=lambda s: s.dropna().astype(str).nunique(),
            target_genes=uniq_join,
        )
        .reset_index()
    )
    genes["n_target_genes"] = genes["target_genes"].map(lambda z: 0 if z == "" else len(z.split(",")))
    tissues = (
        x.groupby(WINDOW_KEY)["tissue"]
        .agg(n_tissues_lead=lambda s: s.dropna().astype(str).nunique())
        .reset_index()
    )
    if "variant_id" in x.columns:
        eqtl_n = (
            x.groupby(WINDOW_KEY)["variant_id"]
            .agg(n_eqtl_lead=lambda s: s.dropna().astype(str).nunique())
            .reset_index()
        )
    else:
        eqtl_n = x.groupby(WINDOW_KEY).size().reset_index(name="n_eqtl_lead")

    out = genes.merge(tissues, on=WINDOW_KEY, how="outer")
    out = out.merge(eqtl_n, on=WINDOW_KEY, how="outer")
    out = out.merge(top, on=WINDOW_KEY, how="left")

    for col in ["n_genes_lead", "n_tissues_lead", "n_eqtl_lead", "n_target_genes"]:
        if col in out.columns:
            out[col] = out[col].fillna(0).astype(int)

    if "target_genes" in out.columns:
        out["target_genes"] = out["target_genes"].fillna("")
    return out

def build_background(windows: pd.DataFrame, lead_eqtl_windows: pd.DataFrame) -> pd.DataFrame:
    bg = windows.copy()
    bg = bg[bg["is_valid"] == 1].copy()
    bg = bg[bg["Fw"] > 0].copy()
    bg = bg.merge(lead_eqtl_windows, on=WINDOW_KEY, how="inner")
    return bg

def build_regions(candidate_windows: pd.DataFrame, max_gap_bp: int) -> pd.DataFrame:
    if candidate_windows.empty:
        return pd.DataFrame(
            columns=[
                "region_id",
                "chrom",
                "region_start",
                "region_end",
                "n_windows",
                "length_bp",
                "Fw_mean",
                "Fw_max",
                "Sw_mean",
                "Sw_max",
                "window_ids",
                "target_genes",
                "n_target_genes",
            ]
        )
    df = candidate_windows.sort_values(["chrom", "win_start", "win_end"], key=chrom_sort_key).reset_index(drop=True)
    chunks = []
    current_idx = [0]

    for i in range(1, len(df)):
        prev = df.iloc[i - 1]
        cur = df.iloc[i]
        same_chrom = cur["chrom"] == prev["chrom"]
        gap = int(cur["win_start"]) - int(prev["win_end"])

        if same_chrom and gap <= max_gap_bp:
            current_idx.append(i)
        else:
            chunks.append(df.iloc[current_idx].copy())
            current_idx = [i]

    chunks.append(df.iloc[current_idx].copy())
    rows = []
    chrom_counter: dict[str, int] = {}

    for chunk in chunks:
        chrom = str(chunk["chrom"].iloc[0])
        chrom_counter[chrom] = chrom_counter.get(chrom, 0) + 1
        region_id = f"{chrom}_region_{chrom_counter[chrom]}"

        region_start = int(chunk["win_start"].min())
        region_end = int(chunk["win_end"].max())
        length_bp = int(region_end - region_start)
        n_windows = int(len(chunk))

        target_genes = union_gene_lists(chunk["target_genes"]) if "target_genes" in chunk.columns else ""
        n_target_genes = 0 if target_genes == "" else len(target_genes.split(","))
        window_ids = ",".join(str(int(x)) for x in chunk["win_id"].tolist())

        rows.append(
            {
                "region_id": region_id,
                "chrom": chunk["chrom"].iloc[0],
                "region_start": region_start,
                "region_end": region_end,
                "n_windows": n_windows,
                "length_bp": length_bp,
                "Fw_mean": float(chunk["Fw"].mean()),
                "Fw_max": float(chunk["Fw"].max()),
                "Sw_mean": float(chunk["Sw_global"].mean()),
                "Sw_max": float(chunk["Sw_global"].max()),
                "window_ids": window_ids,
                "target_genes": target_genes,
                "n_target_genes": n_target_genes,
            }
        )
    return pd.DataFrame(rows)

def write_qc(
    outpath: Path,
    n_input_windows: int,
    fw_q: float,
    sw_q: float,
    fw_thr: float,
    sw_thr: float,
    max_gap_bp: int,
    candidate_windows: pd.DataFrame,
    regions: pd.DataFrame,
) -> None:
    counts = (
        regions.groupby("chrom")
        .size()
        .sort_index(key=lambda x: pd.to_numeric(x, errors="coerce"))
    ) if not regions.empty else pd.Series(dtype=int)

    with outpath.open("w", encoding="utf-8") as f:
        f.write(f"n_input_windows\t{n_input_windows}\n")
        f.write(f"fw_quantile\t{fw_q}\n")
        f.write(f"sw_quantile\t{sw_q}\n")
        f.write(f"Fw_threshold\t{fw_thr}\n")
        f.write(f"Sw_threshold\t{sw_thr}\n")
        f.write(f"max_gap_bp\t{max_gap_bp}\n")
        f.write(f"n_candidate_windows\t{len(candidate_windows)}\n")
        f.write(f"n_candidate_regions\t{len(regions)}\n")
        f.write(f"n_candidate_regions_nw_ge_2\t{int((regions['n_windows'] >= 2).sum()) if not regions.empty else 0}\n")
        f.write(f"n_candidate_regions_nw_ge_5\t{int((regions['n_windows'] >= 5).sum()) if not regions.empty else 0}\n")
        f.write("\n")
        f.write("candidate_region_counts_by_chrom\n")
        f.write("chrom\n")
        if not counts.empty:
            f.write(counts.to_string())
            f.write("\n")

def make_scatter_plot(
    background_windows: pd.DataFrame,
    candidate_windows: pd.DataFrame,
    fw_thr: float,
    sw_thr: float,
    outpath: Path,
) -> None:
    plt.figure(figsize=(8, 6))

    plt.scatter(
        background_windows["Fw"],
        background_windows["Sw_global"],
        s=6,
        alpha=0.20,
        rasterized=True,
        label="Background windows",
    )
    if not candidate_windows.empty:
        plt.scatter(
            candidate_windows["Fw"],
            candidate_windows["Sw_global"],
            s=16,
            alpha=0.85,
            rasterized=True,
            label="Candidate windows",
        )
    plt.axvline(fw_thr, linestyle="--", linewidth=1)
    plt.axhline(sw_thr, linestyle="--", linewidth=1)
    plt.xlabel("Fw")
    plt.ylabel("Sw_global")
    plt.title("Adaptive introgression candidates")
    plt.legend(frameon=False)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300)
    plt.close()

def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    lead_eqtl_path = Path(args.lead_eqtl)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    windows = read_windows(input_path)
    lead_eqtl = read_lead_eqtl(lead_eqtl_path)
    lead_eqtl_windows = lead_eqtl[WINDOW_KEY].drop_duplicates().copy()

    # background = валидный + интрогрессированный + нужный lead-eQTL
    background = build_background(windows, lead_eqtl_windows)
    fw_thr = float(background["Fw"].quantile(args.fw_quantile))
    sw_thr = float(background["Sw_global"].quantile(args.sw_quantile))
    candidate_windows = background[
        (background["Fw"] > fw_thr) & (background["Sw_global"] > sw_thr)
    ].copy()

    # Суммаризация для окон-кандидатов
    candidate_keys = candidate_windows[WINDOW_KEY].drop_duplicates()
    candidate_lead_eqtl = lead_eqtl.merge(candidate_keys, on=WINDOW_KEY, how="inner")
    candidate_lead_summary = summarize_candidate_lead_eqtl(candidate_lead_eqtl)
    candidate_windows = candidate_windows.merge(
        candidate_lead_summary,
        on=WINDOW_KEY,
        how="left",
    )
    for col in ["n_genes_lead", "n_tissues_lead", "n_eqtl_lead", "n_target_genes"]:
        if col in candidate_windows.columns:
            candidate_windows[col] = candidate_windows[col].fillna(0).astype(int)
          
    if "target_genes" in candidate_windows.columns:
        candidate_windows["target_genes"] = candidate_windows["target_genes"].fillna("")
      
    candidate_windows = candidate_windows.sort_values(
        ["chrom", "win_start", "win_end"],
        key=chrom_sort_key
    ).reset_index(drop=True)

    regions = build_regions(candidate_windows, args.max_gap_bp)
    regions_nw_ge_2 = regions[regions["n_windows"] >= 2].copy() if not regions.empty else regions.copy()
    regions_nw_ge_5 = regions[regions["n_windows"] >= 5].copy() if not regions.empty else regions.copy()

    candidate_windows.to_csv(outdir / "adaptive_candidate_windows.tsv", sep="\t", index=False)
    regions.to_csv(outdir / "adaptive_candidate_regions.tsv", sep="\t", index=False)
    regions_nw_ge_2.to_csv(outdir / "adaptive_candidate_regions.nw_ge_2.tsv", sep="\t", index=False)
    regions_nw_ge_5.to_csv(outdir / "adaptive_candidate_regions.nw_ge_5.tsv", sep="\t", index=False)

    write_qc(
        outpath=outdir / "adaptive_candidates_qc.txt",
        n_input_windows=len(background),
        fw_q=args.fw_quantile,
        sw_q=args.sw_quantile,
        fw_thr=fw_thr,
        sw_thr=sw_thr,
        max_gap_bp=args.max_gap_bp,
        candidate_windows=candidate_windows,
        regions=regions,
    )
    make_scatter_plot(
        background_windows=background,
        candidate_windows=candidate_windows,
        fw_thr=fw_thr,
        sw_thr=sw_thr,
        outpath=outdir / "adaptive_introgression_scatter.png",
    )
    print(f"Fw threshold ({args.fw_quantile:.2f}): {fw_thr:.6f}")
    print(f"Sw threshold ({args.sw_quantile:.2f}): {sw_thr:.6f}")
    print(f"Candidate windows: {len(candidate_windows)}")
    print(f"Candidate regions: {len(regions)}")
    print(f"Candidate regions (n_windows >= 2): {len(regions_nw_ge_2)}")
    print(f"Candidate regions (n_windows >= 5): {len(regions_nw_ge_5)}")

    if not regions.empty:
        preview_cols = [
            "region_id",
            "chrom",
            "region_start",
            "region_end",
            "n_windows",
            "Fw_max",
            "Sw_max",
            "target_genes",
        ]
        print(regions[preview_cols].head(20).to_string(index=False))
    return 0

if __name__ == "__main__":
    try:
        raise SystemExit(main())
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        raise SystemExit(1)
