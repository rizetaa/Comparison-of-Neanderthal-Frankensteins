#!/usr/bin/env python3
# Визуализация всех результатов анализа (chr6)
# Пререквизит 03_main_analysis.py

import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
from scipy import stats
import json
import os
import warnings

warnings.filterwarnings("ignore")

WORKDIR = os.path.expanduser("~/nd_pipeline")
OUT_A   = f"{WORKDIR}/results/pipeline_A"
OUT     = f"{WORKDIR}/results/analysis"
FIGS    = f"{WORKDIR}/results/figures"
os.makedirs(FIGS, exist_ok=True)

sns.set_theme(style="whitegrid", font_scale=1.2)
PALETTE = {
    "Zero":       "#CCCCCC",
    "Rare":       "#4575B4",
    "Low":        "#74ADD1",
    "Intermediate":"#FEE090",
    "High":       "#F46D43",
    "Very_High":  "#D73027",
}
FREQ_ORDER = ["Zero", "Rare", "Low", "Intermediate", "High", "Very_High"]
FREQ_ORDER_NONZERO = ["Rare", "Low", "Intermediate", "High", "Very_High"]

print("\n[Fig 1] iSFS — Спектр частот интрогрессии")
isfs = pd.read_csv(f"{OUT}/isfs.tsv", sep="\t")
isfs = isfs.set_index("freq_bin").reindex(FREQ_ORDER).reset_index()
isfs["n_windows"] = isfs["n_windows"].fillna(0).astype(int)

fig, ax = plt.subplots(figsize=(9, 5))
bars = ax.bar(
    isfs["freq_bin"],
    isfs["n_windows"],
    color=[PALETTE.get(b, "#888888") for b in isfs["freq_bin"]],
    edgecolor="black", linewidth=0.7
)
for bar, n in zip(bars, isfs["n_windows"]):
    if n > 0:
        ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 50,
                f"{n:,}", ha="center", va="bottom", fontsize=10)

ax.set_xlabel("Бин частоты интрогрессии (Fw)", fontsize=13)
ax.set_ylabel("Количество окон (1 kb)", fontsize=13)
ax.set_title("Спектр частот неандертальской интрогрессии\n(chr6, окна 1 kb)", fontsize=14)
ax.set_yscale("log")
ax.set_ylim(bottom=1)
plt.tight_layout()
fig.savefig(f"{FIGS}/fig1_iSFS.png", dpi=150, bbox_inches="tight")
fig.savefig(f"{FIGS}/fig1_iSFS.pdf", bbox_inches="tight")
plt.close()
print(f"  Сохранено: {FIGS}/fig1_iSFS.png")

print("[Fig 2] Manhattan plot интрогрессии")
manhattan = pd.read_csv(f"{OUT}/manhattan_data.tsv", sep="\t")
fig, ax = plt.subplots(figsize=(14, 4))
colors_man = [PALETTE.get(b, "#888888") for b in manhattan["freq_bin"]]
ax.scatter(
    manhattan["win_start"] / 1e6,
    manhattan["Fw"],
    c=colors_man, s=1.5, alpha=0.6, linewidths=0
)
ax.set_xlabel("Позиция на chr6 (Mb)", fontsize=13)
ax.set_ylabel("Частота интрогрессии (Fw)", fontsize=13)
ax.set_title("Manhattan plot неандертальской интрогрессии (chr6)", fontsize=14)
ax.set_xlim(0, 171)

legend_patches = [
    mpatches.Patch(color=PALETTE[b], label=b) for b in FREQ_ORDER if b != "Zero"
]
ax.legend(handles=legend_patches, title="Бин частоты",
          loc="upper right", fontsize=9, ncol=2)
plt.tight_layout()
fig.savefig(f"{FIGS}/fig2_manhattan.png", dpi=150, bbox_inches="tight")
fig.savefig(f"{FIGS}/fig2_manhattan.pdf", bbox_inches="tight")
plt.close()
print(f"  Сохранено: {FIGS}/fig2_manhattan.png")

print("[Fig 3] Boxplot Fw_bin vs Sw")
boxplot_data = pd.read_csv(f"{OUT}/boxplot_data.tsv", sep="\t")
boxplot_data = boxplot_data[boxplot_data["freq_bin"].isin(FREQ_ORDER_NONZERO)]
boxplot_data["freq_bin"] = pd.Categorical(
    boxplot_data["freq_bin"], categories=FREQ_ORDER_NONZERO, ordered=True
)
fig, axes = plt.subplots(1, 2, figsize=(14, 6))

# Violin plot
ax = axes[0]
sns.violinplot(
    data=boxplot_data, x="freq_bin", y="Sw_max",
    palette=[PALETTE[b] for b in FREQ_ORDER_NONZERO],
    order=FREQ_ORDER_NONZERO, inner="box", ax=ax, cut=0
)
ax.set_xlabel("Бин частоты интрогрессии", fontsize=12)
ax.set_ylabel("Сила eQTL (Sw = max|Z|)", fontsize=12)
ax.set_title("Violin plot: Fw_bin vs Sw\n(chr6, Fw > 0)", fontsize=13)

medians = boxplot_data.groupby("freq_bin", observed=True)["Sw_max"].median()
for i, b in enumerate(FREQ_ORDER_NONZERO):
    if b in medians.index:
        ax.text(i, medians[b] + 0.1, f"{medians[b]:.2f}",
                ha="center", va="bottom", fontsize=9, color="black", fontweight="bold")
# Boxplot
ax2 = axes[1]
sns.boxplot(
    data=boxplot_data, x="freq_bin", y="Sw_max",
    palette=[PALETTE[b] for b in FREQ_ORDER_NONZERO],
    order=FREQ_ORDER_NONZERO, ax=ax2,
    showfliers=False, width=0.6
)
ax2.set_xlabel("Бин частоты интрогрессии", fontsize=12)
ax2.set_ylabel("Сила eQTL (Sw = max|Z|)", fontsize=12)
ax2.set_title("Boxplot: Fw_bin vs Sw\n(chr6, Fw > 0)", fontsize=13)

x_vals = range(len(FREQ_ORDER_NONZERO))
y_vals = [medians.get(b, np.nan) for b in FREQ_ORDER_NONZERO]
valid = [(x, y) for x, y in zip(x_vals, y_vals) if not np.isnan(y)]
if len(valid) >= 2:
    xv, yv = zip(*valid)
    ax2.plot(xv, yv, "k--o", linewidth=1.5, markersize=5, label="Медиана")
    ax2.legend(fontsize=10)

plt.suptitle("Очищающий отбор: сила eQTL и частота интрогрессии (chr6)",
             fontsize=14, y=1.02)
plt.tight_layout()
fig.savefig(f"{FIGS}/fig3_violin_boxplot.png", dpi=150, bbox_inches="tight")
fig.savefig(f"{FIGS}/fig3_violin_boxplot.pdf", bbox_inches="tight")
plt.close()
print(f"  Сохранено: {FIGS}/fig3_violin_boxplot.png")

print("[Fig 4] Split Violin — Интрогрессия vs Контроль")
violin_data = pd.read_csv(f"{OUT}/violin_data.tsv", sep="\t")
dtss_cats = ["Promoter", "Near", "Distal"]
violin_data = violin_data[violin_data["dtss_cat"].isin(dtss_cats)]

try:
    with open(f"{OUT}/bootstrap_results.json") as f:
        boot_res = json.load(f)
except:
    boot_res = {}

def pval_stars(p):
    if p < 0.001: return "***"
    elif p < 0.01: return "**"
    elif p < 0.05: return "*"
    else: return "ns"

fig, axes = plt.subplots(1, 3, figsize=(15, 6), sharey=True)
colors_split = {"Introgressed": "#E74C3C", "Control": "#3498DB"}

for i, cat in enumerate(dtss_cats):
    ax = axes[i]
    cat_data = violin_data[violin_data["dtss_cat"] == cat]
    
    if len(cat_data) == 0:
        ax.set_title(f"{cat}\n(нет данных)")
        continue
    
    # Split violin
    try:
        # seaborn >= 0.12
        sns.violinplot(
            data=cat_data, x="dtss_cat", y="Sw_max", hue="group",
            inner="box", ax=ax,
            palette=colors_split, cut=0,
            order=[cat], split=True
        )
    except TypeError:
        # seaborn >= 0.13
        sns.violinplot(
            data=cat_data, x="dtss_cat", y="Sw_max", hue="group",
            inner="box", ax=ax,
            palette=colors_split, cut=0,
            order=[cat]
        )
      
    for grp, color in colors_split.items():
        grp_data = cat_data[cat_data["group"] == grp]["Sw_max"].dropna()
        if len(grp_data) > 0:
            med = grp_data.median()
            ax.axhline(med, color=color, linestyle=":", linewidth=1.5, alpha=0.8)
    # p-value из bootstrap
    if cat in boot_res and "pval_bootstrap" in boot_res[cat]:
        p = boot_res[cat]["pval_bootstrap"]
        stars = pval_stars(p)
        n_i = boot_res[cat].get("n_introgressed", "?")
        n_c = boot_res[cat].get("n_control", "?")
        ax.set_title(
            f"{cat}\n{stars} (p={p:.3f})\n"
            f"n_introg={n_i}, n_ctrl={n_c}",
            fontsize=11
        )
    else:
        ax.set_title(f"{cat}", fontsize=11)

    ax.set_xlabel("")
    if i == 0:
        ax.set_ylabel("Сила eQTL (Sw = max|Z|)", fontsize=12)
    else:
        ax.set_ylabel("")
    legend = ax.get_legend()
    if legend is not None:
        legend.remove()

legend_patches = [
    mpatches.Patch(color=colors_split["Introgressed"], label="Интрогрессия (Fw > 0)"),
    mpatches.Patch(color=colors_split["Control"],      label="Контроль (Fw = 0)"),
]
fig.legend(handles=legend_patches, loc="upper center", ncol=2,
           fontsize=11, bbox_to_anchor=(0.5, 1.02))
plt.suptitle("Сравнение силы eQTL: интрогрессия vs контроль\nпо расстоянию до TSS (chr6)",
             fontsize=13, y=1.06)
plt.tight_layout()
fig.savefig(f"{FIGS}/fig4_split_violin.png", dpi=150, bbox_inches="tight")
fig.savefig(f"{FIGS}/fig4_split_violin.pdf", bbox_inches="tight")
plt.close()
print(f"  Сохранено: {FIGS}/fig4_split_violin.png")

print("[Fig 5] Scatter plot Fw vs Sw (адаптивная интрогрессия)")
scatter = pd.read_csv(f"{OUT}/scatter_data.tsv", sep="\t")
try:
    with open(f"{OUT}/thresholds.json") as f:
        thresholds = json.load(f)
    Fw_95 = thresholds["Fw_95"]
    Sw_95 = thresholds["Sw_95"]
except:
    Fw_95 = scatter["Fw"].quantile(0.95)
    Sw_95 = scatter["Sw_max"].quantile(0.95)
fig, ax = plt.subplots(figsize=(10, 7))

non_cand = scatter[scatter["is_candidate"] == 0]
ax.scatter(non_cand["Fw"], non_cand["Sw_max"],
           c="#AAAAAA", s=8, alpha=0.4, linewidths=0, label="Остальные окна")

cand = scatter[scatter["is_candidate"] == 1]
ax.scatter(cand["Fw"], cand["Sw_max"],
           c="#E74C3C", s=40, alpha=0.9, linewidths=0.5,
           edgecolors="darkred", label=f"Кандидаты (n={len(cand)})", zorder=5)

if "nearest_gene" in cand.columns and len(cand) > 0:
    top_cand = cand.nlargest(min(15, len(cand)), "Sw_max")
    for _, row in top_cand.iterrows():
        gene = str(row.get("nearest_gene", ""))
        if gene and gene != "nan":
            ax.annotate(
                gene,
                xy=(row["Fw"], row["Sw_max"]),
                xytext=(5, 3), textcoords="offset points",
                fontsize=7, color="darkred",
                arrowprops=dict(arrowstyle="-", color="gray", lw=0.5)
            )
# Пороговые линии
ax.axvline(Fw_95, color="#E74C3C", linestyle="--", linewidth=1.2,
           label=f"Fw 95% = {Fw_95:.3f}", alpha=0.7)
ax.axhline(Sw_95, color="#3498DB", linestyle="--", linewidth=1.2,
           label=f"Sw 95% = {Sw_95:.2f}", alpha=0.7)

ax.fill_betweenx([Sw_95, scatter["Sw_max"].max() * 1.05],
                  Fw_95, scatter["Fw"].max() * 1.05,
                  alpha=0.08, color="#E74C3C")

ax.set_xlabel("Частота интрогрессии (Fw)", fontsize=13)
ax.set_ylabel("Сила eQTL (Sw = max|Z|)", fontsize=13)
ax.set_title("Адаптивная интрогрессия: Fw vs Sw\n(chr6, топ 5% по обоим критериям)", fontsize=14)
ax.legend(fontsize=10, loc="upper left")
plt.tight_layout()
fig.savefig(f"{FIGS}/fig5_scatter_adaptive.png", dpi=150, bbox_inches="tight")
fig.savefig(f"{FIGS}/fig5_scatter_adaptive.pdf", bbox_inches="tight")
plt.close()
print(f"  Сохранено: {FIGS}/fig5_scatter_adaptive.png")

print("[Fig 6] Сводная панель")
fig = plt.figure(figsize=(16, 12))
gs = fig.add_gridspec(3, 3, hspace=0.45, wspace=0.35)
# iSFS
ax_isfs = fig.add_subplot(gs[0, 0])
isfs_plot = isfs[isfs["freq_bin"] != "Zero"]
ax_isfs.bar(
    isfs_plot["freq_bin"], isfs_plot["n_windows"],
    color=[PALETTE[b] for b in isfs_plot["freq_bin"]],
    edgecolor="black", linewidth=0.5
)
ax_isfs.set_yscale("log")
ax_isfs.set_title("A. iSFS (chr6)", fontsize=11, fontweight="bold")
ax_isfs.set_xlabel("Бин частоты", fontsize=9)
ax_isfs.set_ylabel("Кол-во окон (log)", fontsize=9)
ax_isfs.tick_params(axis="x", rotation=30, labelsize=8)
# 6b: Manhattan
ax_man = fig.add_subplot(gs[0, 1:])
ax_man.scatter(
    manhattan["win_start"] / 1e6, manhattan["Fw"],
    c=[PALETTE.get(b, "#888") for b in manhattan["freq_bin"]],
    s=1, alpha=0.5, linewidths=0
)
ax_man.set_title("B. Manhattan plot (chr6)", fontsize=11, fontweight="bold")
ax_man.set_xlabel("Позиция (Mb)", fontsize=9)
ax_man.set_ylabel("Fw", fontsize=9)
ax_man.set_xlim(0, 171)

# 6c: Boxplot
ax_box = fig.add_subplot(gs[1, :2])
bp_data = boxplot_data[boxplot_data["freq_bin"].isin(FREQ_ORDER_NONZERO)]
sns.boxplot(
    data=bp_data, x="freq_bin", y="Sw_max",
    palette=[PALETTE[b] for b in FREQ_ORDER_NONZERO],
    order=FREQ_ORDER_NONZERO, ax=ax_box, showfliers=False, width=0.6
)
ax_box.set_title("C. Очищающий отбор: Fw_bin vs Sw", fontsize=11, fontweight="bold")
ax_box.set_xlabel("Бин частоты", fontsize=9)
ax_box.set_ylabel("Sw (max|Z|)", fontsize=9)
ax_box.tick_params(axis="x", rotation=20, labelsize=8)

# 6d: Scatter
ax_sc = fig.add_subplot(gs[1, 2])
ax_sc.scatter(non_cand["Fw"], non_cand["Sw_max"],
              c="#CCCCCC", s=3, alpha=0.3, linewidths=0)
ax_sc.scatter(cand["Fw"], cand["Sw_max"],
              c="#E74C3C", s=20, alpha=0.9, linewidths=0, zorder=5)
ax_sc.axvline(Fw_95, color="#E74C3C", linestyle="--", linewidth=1, alpha=0.7)
ax_sc.axhline(Sw_95, color="#3498DB", linestyle="--", linewidth=1, alpha=0.7)
ax_sc.set_title("D. Адаптивная интрогрессия", fontsize=11, fontweight="bold")
ax_sc.set_xlabel("Fw", fontsize=9)
ax_sc.set_ylabel("Sw", fontsize=9)

# 6e: Bootstrap
ax_boot = fig.add_subplot(gs[2, :])
if boot_res:
    cats_boot = [c for c in dtss_cats if c in boot_res and "pval_bootstrap" in boot_res[c]]
    pvals_boot = [boot_res[c]["pval_bootstrap"] for c in cats_boot]
    pvals_mw   = [boot_res[c].get("pval_mw_standard", 1.0) for c in cats_boot]
    x = np.arange(len(cats_boot))
    width = 0.35
    bars1 = ax_boot.bar(x - width/2, [-np.log10(p + 1e-10) for p in pvals_mw],
                         width, label="MW стандартный", color="#3498DB", alpha=0.8)
    bars2 = ax_boot.bar(x + width/2, [-np.log10(p + 1e-10) for p in pvals_boot],
                         width, label="Block Bootstrap", color="#E74C3C", alpha=0.8)
    ax_boot.axhline(-np.log10(0.05), color="black", linestyle="--",
                    linewidth=1, label="p=0.05")
    ax_boot.set_xticks(x)
    ax_boot.set_xticklabels(cats_boot, fontsize=10)
    ax_boot.set_ylabel("-log10(p-value)", fontsize=10)
    ax_boot.set_title("E. Значимость различий Sw: интрогрессия vs контроль (Mann-Whitney)",
                       fontsize=11, fontweight="bold")
    ax_boot.legend(fontsize=9)
    
    for bar, p in zip(bars2, pvals_boot):
        ax_boot.text(bar.get_x() + bar.get_width()/2,
                     bar.get_height() + 0.05,
                     pval_stars(p), ha="center", va="bottom", fontsize=12)

plt.suptitle("Сводная панель: eQTL-влияние неандертальской интрогрессии (chr6)",
             fontsize=14, fontweight="bold", y=1.01)
fig.savefig(f"{FIGS}/fig6_summary_panel.png", dpi=150, bbox_inches="tight")
fig.savefig(f"{FIGS}/fig6_summary_panel.pdf", bbox_inches="tight")
plt.close()
print(f"  Сохранено: {FIGS}/fig6_summary_panel.png")

print("Все рисунки сохранены в:", FIGS)
print("  fig1_iSFS.png/pdf спектр частот интрогрессии")
print("  fig2_manhattan.png/pdf Manhattan plot")
print("  fig3_violin_boxplot.png/pdf Violin/Boxplot Fw_bin, Sw")
print("  fig4_split_violin.png/pdf Split Violin по D_TSS")
print("  fig5_scatter_adaptive.png/pdf Scatter Fw, Sw")
print("  fig6_summary_panel.png/pdf Сводка")
