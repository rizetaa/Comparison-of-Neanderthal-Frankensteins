#!/usr/bin/env python3

from __future__ import annotations
from pathlib import Path
from scipy.stats import spearmanr

import argparse
import sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.formula.api as smf

FREQ_ORDER = ["Rare", "Low", "Intermediate", "High", "Very_High"]

def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Two-part regression analysis for Pipeline A with density control."
    )
    parser.add_argument(
        "--input",
        required=True,
        help="Path to chr*_windows_full.tsv",
    )
    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory",
    )
    parser.add_argument(
        "--plot-type",
        default="box",
        choices=["box", "violin"],
        help="Plot type for Sw by frequency bin",
    )
    return parser.parse_args()

def load_data(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t")
    required = {
        "win_id",
        "Fw",
        "Sw_global",
        "DTSS",
        "recomb_rate_cM_Mb",
        "freq_bin",
        "seg_id",
        "n_eqtl_lead",
    }
    missing = required - set(df.columns)
    if missing:
        raise RuntimeError(f"Missing required columns in {path}: {sorted(missing)}")
    return df

def prepare_tables(df: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    # Только интрогрессированные окна
    df = df[df["Fw"] > 0].copy()
    if df.empty:
        raise RuntimeError("No windows with Fw > 0")

    # Бинарный вывод
    df["Sw_binary"] = (df["Sw_global"] > 0).astype(int)
    df["log_DTSS"] = np.log1p(df["DTSS"])
    df["log_recomb"] = np.log1p(df["recomb_rate_cM_Mb"])
    df["log_n_eqtl"] = np.log1p(df["n_eqtl_lead"])
    # Только ячейки для анализа
    df = df[df["freq_bin"].isin(FREQ_ORDER)].copy()
    df["freq_bin"] = pd.Categorical(df["freq_bin"], categories=FREQ_ORDER, ordered=True)
    #  Идентификаторы кластера для кластеризованного SE
    df = df[df["seg_id"].notna()].copy()
    if df.empty:
        raise RuntimeError("No windows with non-missing seg_id after filtering")

    # Шаг 2: только окна с сигналом
    df_lin = df[df["Sw_global"] > 0].copy()
    if df_lin.empty:
        raise RuntimeError("No windows with Sw_global > 0 among Fw > 0 windows")
    return df, df_lin

def fit_logistic(df: pd.DataFrame):
    model = smf.logit(
        "Sw_binary ~ Fw + log_DTSS + log_recomb",
        data=df,
    ).fit(
        disp=False,
        cov_type="cluster",
        cov_kwds={"groups": df["seg_id"]},
    )
    return model

def fit_linear(df_lin: pd.DataFrame):
    model = smf.ols(
        "Sw_global ~ Fw + log_DTSS + log_recomb + log_n_eqtl",
        data=df_lin,
    ).fit(
        cov_type="cluster",
        cov_kwds={"groups": df_lin["seg_id"]},
    )
    return model

def save_summary(model, path: Path, title: str) -> None:
    with path.open("w", encoding="utf-8") as f:
        f.write(title + "\n")
        f.write("=" * len(title) + "\n\n")
        f.write(model.summary().as_text())
        f.write("\n")

def summarize_medians(df_lin: pd.DataFrame) -> pd.DataFrame:
    med = (
        df_lin.groupby("freq_bin", observed=False)["Sw_global"]
        .median()
        .reindex(FREQ_ORDER)
        .reset_index()
        .rename(columns={"Sw_global": "median_Sw"})
    )
    med["bin_index"] = np.arange(1, len(med) + 1)
    return med

def make_plot(df_lin: pd.DataFrame, out_path: Path, plot_type: str) -> tuple[float, float]:
    plot_df = df_lin.copy()
    plot_df = plot_df.sort_values("freq_bin")
    grouped = [
        plot_df.loc[plot_df["freq_bin"] == b, "Sw_global"].values
        for b in FREQ_ORDER
    ]
    plt.figure(figsize=(8, 5))
    x = np.arange(1, len(FREQ_ORDER) + 1)

    if plot_type == "violin":
        plt.violinplot(grouped, positions=x, showmedians=True)
        plt.xticks(x, FREQ_ORDER, rotation=20)
    else:
        plt.boxplot(grouped, tick_labels=FREQ_ORDER)
        plt.xticks(rotation=20)

    med = summarize_medians(plot_df)
    plt.plot(x, med["median_Sw"].values, marker="o")
    rho, pval = spearmanr(med["bin_index"], med["median_Sw"])
    plt.xlabel("Frequency bin")
    plt.ylabel("Sw")
    plt.title(f"Sw by introgression frequency bin\nSpearman rho={rho:.3f}, p={pval:.3g}")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()
    return float(rho), float(pval)

def save_qc(df: pd.DataFrame, df_lin: pd.DataFrame, out_path: Path) -> None:
    lines = []
    lines.append(f"n_windows_Fw_gt_0\t{len(df)}")
    lines.append(f"n_windows_Sw_gt_0\t{len(df_lin)}")
    lines.append(f"n_clusters\t{df['seg_id'].nunique()}")
    lines.append("")
    lines.append("freq_bin_counts")
    lines.append(df["freq_bin"].value_counts().reindex(FREQ_ORDER).to_string())
    lines.append("")
    lines.append("Sw_summary_positive_only")
    lines.append(df_lin["Sw_global"].describe().to_string())
    lines.append("")
    lines.append("n_eqtl_lead_summary")
    lines.append(df["n_eqtl_lead"].describe().to_string())
    out_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

def save_median_table(df_lin: pd.DataFrame, out_path: Path) -> tuple[float, float]:
    med = summarize_medians(df_lin)
    rho, pval = spearmanr(med["bin_index"], med["median_Sw"])
    med["spearman_rho"] = rho
    med["spearman_pvalue"] = pval
    med.to_csv(out_path, sep="\t", index=False)
    return float(rho), float(pval)

def main() -> int:
    args = parse_args()
    input_path = Path(args.input)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = load_data(input_path)
    df, df_lin = prepare_tables(df)
    # Сохраняем таблицы
    df.to_csv(outdir / "analysis_windows_Fw_gt_0.tsv", sep="\t", index=False)
    df_lin.to_csv(outdir / "analysis_windows_Fw_gt_0_Sw_gt_0.tsv", sep="\t", index=False)
    # Модели
    logit_model = fit_logistic(df)
    linear_model = fit_linear(df_lin)
    save_summary(
        logit_model,
        outdir / "logistic_model_summary.txt",
        "Step 1: Logistic regression (Sw > 0) with density control",
    )
    save_summary(
        linear_model,
        outdir / "linear_model_summary.txt",
        "Step 2: Linear regression (Sw | Sw > 0) with density control",
    )
    # Таблица коэффициентов
    coef_rows = []
    for model_name, model in [("logistic", logit_model), ("linear", linear_model)]:
        ci = model.conf_int()
        for term in model.params.index:
            coef_rows.append({
                "model": model_name,
                "term": term,
                "coef": model.params[term],
                "se": model.bse[term],
                "pvalue": model.pvalues[term],
                "ci_low": ci.loc[term, 0],
                "ci_high": ci.loc[term, 1],
            })
    pd.DataFrame(coef_rows).to_csv(outdir / "model_coefficients.tsv", sep="\t", index=False)
  
    # Медианы + график
    rho, pval = save_median_table(df_lin, outdir / "Sw_medians_by_freq_bin.tsv")
    make_plot(df_lin, outdir / f"Sw_by_freq_bin_{args.plot_type}.png", args.plot_type)
    # qc
    save_qc(df, df_lin, outdir / "analysis_qc.txt")
    print("Logistic Fw coef:", round(float(logit_model.params["Fw"]), 6))
    print("Logistic Fw pvalue:", f"{float(logit_model.pvalues['Fw']):.3g}")
    print("Linear Fw coef:", round(float(linear_model.params["Fw"]), 6))
    print("Linear Fw pvalue:", f"{float(linear_model.pvalues['Fw']):.3g}")
    print("Spearman rho:", round(rho, 6))
    print("Spearman pvalue:", f"{pval:.3g}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
