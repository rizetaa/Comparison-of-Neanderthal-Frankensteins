#!/bin/bash

set -euo pipefail
WORKDIR="$HOME/nd_pipeline"
DATA="$WORKDIR/data"
RES="$WORKDIR/results"
LOG="$WORKDIR/logs/pipeline_A.log"

mkdir -p "$RES/pipeline_A"
exec > >(tee -a "$LOG") 2>&1
echo "[$(date)] Pipeline A, старт"

# Параметры
CHR=6
WIN_SIZE=1000
NIS="$DATA/raw/IBS.YRI.grch37.chr6.em.tsv"
HAPMAP="$DATA/hapmap/hapmap_chr6.txt"
GENCODE_GTF="$DATA/gencode/gencode.v19.annotation.gtf"
GTEX_DIR="$DATA/gtex/GTEx_Analysis_v7_eQTL"
OUT="$RES/pipeline_A"
CHR6_LEN=171115067

# Карта рекомбинации HapMap
# Векторный подход через numpy
echo "[$(date)] Шаг 8: Вычисление средней скорости рекомбинации для окон"

python3 << 'PYEOF'
import pandas as pd
import numpy as np
import os

WORKDIR = os.path.expanduser("~/nd_pipeline")
OUT = f"{WORKDIR}/results/pipeline_A"
HAPMAP = f"{WORKDIR}/data/hapmap/hapmap_chr6.txt"

# Читаем HapMap
hapmap = pd.read_csv(HAPMAP, sep=r"\s+", comment="#")
hapmap.columns = [c.strip() for c in hapmap.columns]

# Нормализуем имена колонок
col_map = {}
for c in hapmap.columns:
    cl = c.lower()
    if "position" in cl or "pos" in cl:
        col_map[c] = "pos"
    elif "rate" in cl:
        col_map[c] = "rate"
    elif "map" in cl:
        col_map[c] = "map_cm"
        
hapmap = hapmap.rename(columns=col_map)
hapmap = hapmap.sort_values("pos").reset_index(drop=True)
print(f"  HapMap записей: {len(hapmap)}")
print(f"  Колонки: {list(hapmap.columns)}")
windows = pd.read_csv(f"{OUT}/chr6_windows_full.tsv", sep="\t")

# Векторный расчет средней скорости рекомбинации
# Для каждого окна [win_start, win_end) HapMap точки внутри
# Средний rate, если точек нет, интерполяция из ближайших
pos_arr  = hapmap["pos"].values
rate_arr = hapmap["rate"].values
win_starts = windows["win_start"].values
win_ends   = windows["win_end"].values

# Индексы левой и правой границы для каждого окна
left_idx  = np.searchsorted(pos_arr, win_starts, side="left")
right_idx = np.searchsorted(pos_arr, win_ends,   side="left")
recomb_rates = np.empty(len(windows), dtype=float)

for i in range(len(windows)):
    l, r = left_idx[i], right_idx[i]
    if l < r:
        # Есть HapMap точки внутри окна, значит среднее
        recomb_rates[i] = rate_arr[l:r].mean()
    else:
        if l == 0:
            recomb_rates[i] = rate_arr[0]
        elif l >= len(pos_arr):
            recomb_rates[i] = rate_arr[-1]
        else:
            recomb_rates[i] = (rate_arr[l - 1] + rate_arr[l]) / 2.0

# Замена нулей и отрицательных значений на мининимальное положительное
min_rate = recomb_rates[recomb_rates > 0].min() if (recomb_rates > 0).any() else 1e-4
recomb_rates = np.where(recomb_rates <= 0, min_rate, recomb_rates)
windows["recomb_rate"] = recomb_rates

windows.to_csv(f"{OUT}/chr6_windows_full.tsv", sep="\t", index=False)
print(f"  Обновлена матрица с recomb_rate: {OUT}/chr6_windows_full.tsv")
print(f"  Средняя скорость рекомбинации: {windows['recomb_rate'].mean():.4f} cM/Mb")
PYEOF

echo "[$(date)] Шаг 8 завершен"
echo "[$(date)] Pipeline A завершен"
echo "  Итоговый файл: $OUT/chr6_windows_full.tsv"
echo "  Колонки: chr, win_start, win_end, win_id, Fw, n_hap_introgressed,"
echo "           freq_bin, Sw_max, has_eqtl, D_TSS, dtss_cat, nearest_gene,"
echo "           gene_id, recomb_rate"
