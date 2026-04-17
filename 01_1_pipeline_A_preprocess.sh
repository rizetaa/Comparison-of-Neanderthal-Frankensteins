#!/bin/bash
# Pipeline A, препроцессинг для основного анализа
# Пререквизиты: bedtools, python3 (pandas, numpy), tabix
# Запуск через run

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

# Вычисление Fw
# Нарезка генома chr6 на окна 1000 bp
echo "[$(date)] Шаг 1: Создание окон 1000 bp для chr6"

# Для bedtools
echo -e "6\t${CHR6_LEN}" > "$OUT/chr6.genome"
bedtools makewindows \
    -g "$OUT/chr6.genome" \
    -w $WIN_SIZE \
    > "$OUT/chr6_windows_1kb.bed"
echo "  Создано окон $(wc -l < "$OUT/chr6_windows_1kb.bed")"

# Конвертируем NIS в BED и вычисляем Fw
echo "[$(date)] Шаг 2: Вычисление глобальной частоты интрогрессии Fw"
# Конвертируем NIS в BED и считаем H
python3 << 'PYEOF'
import pandas as pd
import numpy as np
import os

WORKDIR = os.path.expanduser("~/nd_pipeline")
OUT = f"{WORKDIR}/results/pipeline_A"

# Sample CHROM Start End Length
# Start и End, 1-based, включительно
nis = pd.read_csv(f"{WORKDIR}/data/raw/IBS.YRI.grch37.chr6.em.tsv", sep=r"\s+")
print(f"  NIS записей: {len(nis)}")
print(f"  Колонки: {list(nis.columns)}")
print(f"  Пример:\n{nis.head(2).to_string()}")

nis["sample_id"] = nis["Sample"].str.rsplit("_", n=1).str[0]
nis["hap"] = nis["Sample"].str.rsplit("_", n=1).str[-1]
# Уникальные гаплотипы
all_haplotypes = nis["Sample"].unique()
H = len(all_haplotypes)
print(f"  Всего гаплотипов H = {H}")
# Сохраняем H для следующих шагов
with open(f"{OUT}/H_total.txt", "w") as f:
    f.write(str(H))

# Конвертируем NIS в BED (0-based half-open)
# Start и End, 1-based включительно и start = Start - 1, end = End
nis_bed = nis[["Start", "End", "Sample"]].copy()
nis_bed.insert(0, "chr", "6")
nis_bed["bed_start"] = nis_bed["Start"] - 1 # 1-based -> 0-based
nis_bed["bed_end"]   = nis_bed["End"] # включительный -> не включительный (End = last_pos)
# Формат BED: chr bed_start bed_end haplotype_id
nis_bed[["chr", "bed_start", "bed_end", "Sample"]].to_csv(
    f"{OUT}/nis_segments.bed", sep="\t", index=False, header=False)
print(f"  NIS BED сохранен: {OUT}/nis_segments.bed")

# Медианная длина NIS-сегмента
median_seg_len = nis["Length"].median()
print(f"  Медианная длина NIS-сегмента: {median_seg_len:.0f} bp")
with open(f"{OUT}/median_seg_len.txt", "w") as f:
    f.write(str(int(median_seg_len)))
PYEOF

# 2b: bedtools intersect: окна x NIS-сегменты
# фильтруем >= 500 bp
echo "  bedtools intersect: окна × NIS-сегменты"
sort -k1,1 -k2,2n "$OUT/chr6_windows_1kb.bed"  > "$OUT/chr6_windows_sorted.bed"
sort -k1,1 -k2,2n "$OUT/nis_segments.bed"       > "$OUT/nis_segments_sorted.bed"

bedtools intersect \
    -a "$OUT/chr6_windows_sorted.bed" \
    -b "$OUT/nis_segments_sorted.bed" \
    -wo \
    > "$OUT/windows_x_nis.bed"
echo "  Пересечений окна x NIS: $(wc -l < "$OUT/windows_x_nis.bed")"

#2c: вычисляем Fw векторно
python3 << 'PYEOF'
import pandas as pd
import numpy as np
import os

WORKDIR = os.path.expanduser("~/nd_pipeline")
OUT = f"{WORKDIR}/results/pipeline_A"
WIN_SIZE = 1000
MIN_OVERLAP = WIN_SIZE // 2  # 500 bp

# Читаем H
with open(f"{OUT}/H_total.txt") as f:
    H = int(f.read().strip())

# Читаем окна из chrN_windows_sorted.bed (chr, start, end)
# win_id = порядковый номер строки
windows = pd.read_csv(f"{OUT}/chr6_windows_sorted.bed", sep="\t",
                      header=None, names=["chr","win_start","win_end"])
windows["win_id"] = range(len(windows))
windows["chr"] = windows["chr"].astype(str)

# Читаем пересечения
# bedtools intersect 3 колонки + 4 колонки NIS + overlap_bp = 8 колонок
# Окна из chr6_windows_sorted.bed (chr,start,end)
# nis_segments_sorted.bed (chr,start,end,hap_id)
cols = ["chr_win","win_start","win_end",
        "chr_nis","seg_start","seg_end","hap_id","overlap_bp"]
ix = pd.read_csv(f"{OUT}/windows_x_nis.bed", sep="\t",
                 header=None, names=cols)
ix["chr_win"] = ix["chr_win"].astype(str)
print(f"  Пересечений загружено: {len(ix)}")

# Перекрытие >= 500 bp
ix_filt = ix[ix["overlap_bp"] >= MIN_OVERLAP].copy()
print(f"  После фильтра (>=500 bp): {len(ix_filt)}")

# Присваиваем win_id через merge по координатам
ix_filt = ix_filt.merge(
    windows[["chr","win_start","win_end","win_id"]].rename(columns={"chr":"chr_win"}),
    on=["chr_win","win_start","win_end"],
    how="left"
)
# Для каждого окна считаем уникальные гаплотипы с перекрытием >= 500 bp
n_introgressed = (
    ix_filt.groupby("win_id")["hap_id"]
    .nunique()
    .reset_index()
    .rename(columns={"hap_id": "n_hap_introgressed"})
)
# Объединяем с полным списком окон (Fw=0 для окон без пересечений)
fw_df = windows.merge(n_introgressed, on="win_id", how="left")
fw_df["n_hap_introgressed"] = fw_df["n_hap_introgressed"].fillna(0).astype(int)
fw_df["Fw"] = fw_df["n_hap_introgressed"] / H

# Присваиваем бины частот
def assign_freq_bin(fw):
    if fw == 0:
        return "Zero"
    elif fw <= 0.02:
        return "Rare"
    elif fw <= 0.05:
        return "Low"
    elif fw <= 0.10:
        return "Intermediate"
    elif fw <= 0.20:
        return "High"
    else:
        return "Very_High"
        
fw_df["freq_bin"] = fw_df["Fw"].apply(assign_freq_bin)

# Убеждаемся что chr есть
if "chr" not in fw_df.columns:
    fw_df.insert(0, "chr", 6)
else:
    fw_df["chr"] = 6

fw_df.to_csv(f"{OUT}/chr6_windows_Fw.tsv", sep="\t", index=False)
print(f"  Сохранено: {OUT}/chr6_windows_Fw.tsv")
print(f"  Распределение по бинам:\n{fw_df['freq_bin'].value_counts()}")
# BED для bedtools intersect
fw_df[["chr","win_start","win_end","win_id"]].to_csv(
    f"{OUT}/chr6_windows.bed", sep="\t", index=False, header=False)
PYEOF

echo "[$(date)] Шаг 2 завершен"
