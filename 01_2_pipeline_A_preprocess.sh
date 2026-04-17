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

# eQTL
# Обработка GTEx eQTL, фильтрация и пересечение с окнами
echo "[$(date)] Шаг 3: Обработка GTEx eQTL"

python3 << 'PYEOF'
import pandas as pd
import numpy as np
import os
import glob

WORKDIR = os.path.expanduser("~/nd_pipeline")
GTEX_DIR = f"{WORKDIR}/data/gtex/GTEx_Analysis_v7_eQTL"
OUT = f"{WORKDIR}/results/pipeline_A"

# Читаем значимые пары GTEx
eqtl_files = glob.glob(f"{GTEX_DIR}/*.signif_variant_gene_pairs.txt.gz")
if not eqtl_files:
    eqtl_files = glob.glob(f"{GTEX_DIR}/*.signif_variant_gene_pairs.txt")
print(f"  Найдено файлов GTEx: {len(eqtl_files)}")

all_eqtl = []
for f in eqtl_files:
    tissue = os.path.basename(f).replace(".v7.signif_variant_gene_pairs.txt.gz","") \
                                .replace(".v7.signif_variant_gene_pairs.txt","")
    try:
        df = pd.read_csv(f, sep="\t", compression="gzip" if f.endswith(".gz") else None)
    except Exception as e:
        print(f"  Ошибка чтения {f}: {e}")
        continue
    df["tissue"] = tissue
    all_eqtl.append(df)

if not all_eqtl:
    print("  GTEx файлы не найдены, проверьте путь")
    exit(1)

eqtl = pd.concat(all_eqtl, ignore_index=True)
print(f"  Всего eQTL записей: {len(eqtl)}")
# Смотрим формат variant_id и уникальные хромосомы
print(f"  Примеры variant_id:\n{eqtl['variant_id'].head(5).to_string()}")
chrom_samples = eqtl["variant_id"].str.split("_").str[0].value_counts().head(10)
print(f"  Уникальные хромосомы в variant_id (10):\n{chrom_samples.to_string()}")

# Парсим variant_id
# chr6_1000_A_G_b37 иногда без "chr"
# Разбиваем по "_" с maxsplit=4
split_df = eqtl["variant_id"].str.split("_", n=4, expand=True)
split_df.columns = ["vchr","vpos","vref","valt","vbuild"]
eqtl = pd.concat([eqtl, split_df], axis=1)

# Фильтр chr6 или 6
eqtl = eqtl[eqtl["vchr"].isin(["chr6", "6"])].copy()
print(f"  eQTL на chr6: {len(eqtl)}")
print(f"  Уникальные значения vchr: {eqtl['vchr'].unique()[:10] if len(eqtl) > 0 else 'нет данных'}")

# Если 0, то фильтр по строке
if len(eqtl) == 0:
    eqtl_all = pd.concat(all_eqtl, ignore_index=True)
    # Фильтруем по variant_id
    mask = eqtl_all["variant_id"].str.startswith("chr6_") | \
           eqtl_all["variant_id"].str.startswith("6_")
    eqtl = eqtl_all[mask].copy()
    print(f"  После фильтра по строке: {len(eqtl)}")
    if len(eqtl) > 0:
        split_df2 = eqtl["variant_id"].str.split("_", n=4, expand=True)
        split_df2.columns = ["vchr","vpos","vref","valt","vbuild"]
        eqtl = pd.concat([eqtl.drop(columns=["vchr","vpos","vref","valt","vbuild"],
                                     errors="ignore"), split_df2], axis=1)

# Фильтр биаллельные SNP (ref и alt - одна буква)
eqtl = eqtl[
    (eqtl["vref"].str.len() == 1) &
    (eqtl["valt"].str.len() == 1)
].copy()
print(f"  После фильтра биаллельных: {len(eqtl)}")
eqtl["vpos"] = eqtl["vpos"].astype(int)

# Z-score
eqtl["Z"] = eqtl["slope"] / eqtl["slope_se"]
eqtl["absZ"] = eqtl["Z"].abs()

eqtl.to_csv(f"{OUT}/chr6_eqtl_all.tsv", sep="\t", index=False)
print(f"  Сохранен полный список: {OUT}/chr6_eqtl_all.tsv")

if len(eqtl) == 0:
    print("  Предупреждение: eQTL chr6 не найдены, создаем пустой BED.")
    # Пустой файл чтобы bedtools не упал
    open(f"{OUT}/chr6_eqtl.bed", "w").close()
else:
    # BED для пересечения с окнами (0-based: pos-1, pos)
    eqtl_bed = eqtl[["vchr","vpos","gene_id","tissue","slope","slope_se","Z","absZ","variant_id"]].copy()
    eqtl_bed["vstart"] = eqtl_bed["vpos"] - 1
    eqtl_bed["vend"]   = eqtl_bed["vpos"]
    eqtl_bed["vchr_num"] = "6"  # без chr для совместимости с NIS
    eqtl_bed = eqtl_bed[["vchr_num","vstart","vend","variant_id","gene_id","tissue","slope","slope_se","Z","absZ"]]
    eqtl_bed.to_csv(f"{OUT}/chr6_eqtl.bed", sep="\t", index=False, header=False)
    print(f"  BED файл eQTL: {OUT}/chr6_eqtl.bed ({len(eqtl_bed)} записей)")
PYEOF

echo "[$(date)] Шаг 3 завершен"

# bedtools intersect: eQTL SNP x окна
echo "[$(date)] Шаг 4: Пересечение eQTL с окнами (bedtools intersect)"
sort -k1,1 -k2,2n "$OUT/chr6_windows.bed" > "$OUT/chr6_windows_sorted2.bed"
sort -k1,1 -k2,2n "$OUT/chr6_eqtl.bed"   > "$OUT/chr6_eqtl_sorted.bed"

# Пересечение для каждого eQTL SNP находим окно
bedtools intersect \
    -a "$OUT/chr6_eqtl_sorted.bed" \
    -b "$OUT/chr6_windows_sorted2.bed" \
    -wa -wb \
    > "$OUT/chr6_eqtl_x_windows.bed"

echo "  Пересечений: $(wc -l < "$OUT/chr6_eqtl_x_windows.bed")"
