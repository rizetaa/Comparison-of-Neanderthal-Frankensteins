#!/bin/bash
# Pipeline A, препроцессинг для основного анализа
# Надо установить bedtools, python3 (pandas, numpy), tabix
# Запуск через SLURM - sbatch run_01.sh

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

# Нарезка генома chr6 на окна 1000 bp
echo "[$(date)] Шаг 1: Создание окон 1000 bp для chr6"

# Создаём genome file для bedtools
echo -e "6\t${CHR6_LEN}" > "$OUT/chr6.genome"
bedtools makewindows \
    -g "$OUT/chr6.genome" \
    -w $WIN_SIZE \
    > "$OUT/chr6_windows_1kb.bed"
echo "  Создано окон $(wc -l < "$OUT/chr6_windows_1kb.bed")"

# Конвертация NIS в BED и вычисление глобальной частоты Fw
echo "[$(date)] Шаг 2: Вычисление глобальной частоты интрогрессии Fw"
# Конвертируем NIS в BED и считаем H
python3 << 'PYEOF'
import pandas as pd
import numpy as np
import os

WORKDIR = os.path.expanduser("~/nd_pipeline")
OUT = f"{WORKDIR}/results/pipeline_A"

# Читаем tsv-файл
# Sample CHROM Start End Length
# Start и End — 1-based, включительно
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
# Start и End — 1-based включительно и start = Start - 1 (0-based), end = End (не включительно)
nis_bed = nis[["Start", "End", "Sample"]].copy()
nis_bed.insert(0, "chr", "6")
nis_bed["bed_start"] = nis_bed["Start"] - 1   # 1-based → 0-based
nis_bed["bed_end"]   = nis_bed["End"]          # включительный -> не включительный (End уже = last_pos)
# Формат BED: chr bed_start bed_end haplotype_id
nis_bed[["chr", "bed_start", "bed_end", "Sample"]].to_csv(
    f"{OUT}/nis_segments.bed", sep="\t", index=False, header=False)
print(f"  NIS BED сохранен: {OUT}/nis_segments.bed")

# Медианная длина NIS-сегмента (для bootstrap)
median_seg_len = nis["Length"].median()
print(f"  Медианная длина NIS-сегмента: {median_seg_len:.0f} bp")
with open(f"{OUT}/median_seg_len.txt", "w") as f:
    f.write(str(int(median_seg_len)))
PYEOF

# 2b: bedtools intersect: окна x NIS-сегменты
# wo выводит перекрытие в bp, фильтруем >= 500 bp (50% от 1000 bp)
echo "  bedtools intersect: окна × NIS-сегменты"
sort -k1,1 -k2,2n "$OUT/chr6_windows_1kb.bed"  > "$OUT/chr6_windows_sorted.bed"
sort -k1,1 -k2,2n "$OUT/nis_segments.bed"       > "$OUT/nis_segments_sorted.bed"

bedtools intersect \
    -a "$OUT/chr6_windows_sorted.bed" \
    -b "$OUT/nis_segments_sorted.bed" \
    -wo \
    > "$OUT/windows_x_nis.bed"
echo "  Пересечений окна x NIS: $(wc -l < "$OUT/windows_x_nis.bed")"

#2c: Вычисляем Fw векторно
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

# Читаем окна из chr6_windows_sorted.bed (chr, start, end)
# win_id = порядковый номер строки
windows = pd.read_csv(f"{OUT}/chr6_windows_sorted.bed", sep="\t",
                      header=None, names=["chr","win_start","win_end"])
windows["win_id"] = range(len(windows))
windows["chr"] = windows["chr"].astype(str)

# Читаем пересечения
# bedtools intersect 3 колонки окна + 4 колонки NIS + overlap_bp = 8 колонок
# Окна из chr6_windows_sorted.bed (chr,start,end)
# NIS из nis_segments_sorted.bed (chr,start,end,hap_id)
cols = ["chr_win","win_start","win_end",
        "chr_nis","seg_start","seg_end","hap_id","overlap_bp"]
ix = pd.read_csv(f"{OUT}/windows_x_nis.bed", sep="\t",
                 header=None, names=cols)
ix["chr_win"] = ix["chr_win"].astype(str)
print(f"  Пересечений загружено: {len(ix)}")

# Фильтр: перекрытие >= 500 bp
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

# Убеждаемся что chr есть (нормализуем)
if "chr" not in fw_df.columns:
    fw_df.insert(0, "chr", 6)
else:
    fw_df["chr"] = 6

# Сохраняем
fw_df.to_csv(f"{OUT}/chr6_windows_Fw.tsv", sep="\t", index=False)
print(f"  Сохранено: {OUT}/chr6_windows_Fw.tsv")
print(f"  Распределение по бинам:\n{fw_df['freq_bin'].value_counts()}")
# Сохраняем BED для bedtools intersect (только числовой хром)
fw_df[["chr","win_start","win_end","win_id"]].to_csv(
    f"{OUT}/chr6_windows.bed", sep="\t", index=False, header=False)
PYEOF
echo "[$(date)] Шаг 2 завершен"

# Обработка GTEx eQTL — фильтрация и пересечение с окнами
echo "[$(date)] Шаг 3: Обработка GTEx eQTL"

python3 << 'PYEOF'
import pandas as pd
import numpy as np
import os
import glob

WORKDIR = os.path.expanduser("~/nd_pipeline")
GTEX_DIR = f"{WORKDIR}/data/gtex/GTEx_Analysis_v7_eQTL"
OUT = f"{WORKDIR}/results/pipeline_A"
# Читаем значимые пары GTEx v7
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
# Диагностика: смотрим реальный формат variant_id и уникальные хромосомы
print(f"  Примеры variant_id:\n{eqtl['variant_id'].head(5).to_string()}")
# Показываем уникальные первые части из variant_id
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
print(f"  Уникальные значения vchr (10): {eqtl['vchr'].unique()[:10] if len(eqtl) > 0 else 'нет данных'}")

# Если 0 — фильтр по строке
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

# Фильтр биаллельные SNP (ref и alt — одна буква)
eqtl = eqtl[
    (eqtl["vref"].str.len() == 1) &
    (eqtl["valt"].str.len() == 1)
].copy()
print(f"  После фильтра биаллельных: {len(eqtl)}")

eqtl["vpos"] = eqtl["vpos"].astype(int)

# Z-score
eqtl["Z"] = eqtl["slope"] / eqtl["slope_se"]
eqtl["absZ"] = eqtl["Z"].abs()

# Сохраняем полный список eQTL chr6
eqtl.to_csv(f"{OUT}/chr6_eqtl_all.tsv", sep="\t", index=False)
print(f"  Сохранен полный список: {OUT}/chr6_eqtl_all.tsv")

if len(eqtl) == 0:
    print("  warning: eQTL chr6 не найдены — создаём пустой BED.")
    # Пустой файл чтобы bedtools не упал
    open(f"{OUT}/chr6_eqtl.bed", "w").close()
else:
    # BED для пересечения с окнами (0-based: pos-1, pos)
    eqtl_bed = eqtl[["vchr","vpos","gene_id","tissue","slope","slope_se","Z","absZ","variant_id"]].copy()
    eqtl_bed["vstart"] = eqtl_bed["vpos"] - 1
    eqtl_bed["vend"]   = eqtl_bed["vpos"]
    eqtl_bed["vchr_num"] = "6"  # без "chr" для совместимости с NIS
    eqtl_bed = eqtl_bed[["vchr_num","vstart","vend","variant_id","gene_id","tissue","slope","slope_se","Z","absZ"]]
    eqtl_bed.to_csv(f"{OUT}/chr6_eqtl.bed", sep="\t", index=False, header=False)
    print(f"  BED файл eQTL: {OUT}/chr6_eqtl.bed ({len(eqtl_bed)} записей)")
PYEOF
echo "[$(date)] Шаг 3 завершен"

# bedtools intersect — eQTL SNP x окна
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

# Proxy-clumping и агрегация Sw,t
# Seg_assignments через bedtools
echo "[$(date)] Шаг 5: Proxy-clumping и вычисление Sw,t"

python3 << 'PYEOF'
import pandas as pd
import numpy as np
import os

WORKDIR = os.path.expanduser("~/nd_pipeline")
OUT = f"{WORKDIR}/results/pipeline_A"

# Результат пересечения
cols = ["chr_eqtl","start_eqtl","end_eqtl","variant_id","gene_id","tissue",
        "slope","slope_se","Z","absZ",
        "chr_win","win_start","win_end","win_id"]
df = pd.read_csv(f"{OUT}/chr6_eqtl_x_windows.bed", sep="\t",
                 header=None, names=cols)
print(f"  Записей после пересечения: {len(df)}")

# Proxy-clumping (win_id, tissue, gene_id)
# SNP с максимальным |Z|
df["absZ"] = df["Z"].abs()
lead = df.loc[df.groupby(["win_id","tissue","gene_id"])["absZ"].idxmax()].copy()
print(f"  Lead SNP после clumping: {len(lead)}")

# Агрегация Sw,t = max_g |Z_lead|
sw = lead.groupby(["win_id","tissue"])["absZ"].max().reset_index()
sw.columns = ["win_id","tissue","Sw"]

# Sw = max по всем тканям для каждого окна
sw_max = sw.groupby("win_id")["Sw"].max().reset_index()
sw_max.columns = ["win_id","Sw_max"]

# Объединяем с Fw
fw_df = pd.read_csv(f"{OUT}/chr6_windows_Fw.tsv", sep="\t")
fw_sw = fw_df.merge(sw_max, on="win_id", how="left")
fw_sw["Sw_max"] = fw_sw["Sw_max"].fillna(0)
fw_sw["has_eqtl"] = (fw_sw["Sw_max"] > 0).astype(int)

fw_sw.to_csv(f"{OUT}/chr6_windows_Fw_Sw.tsv", sep="\t", index=False)
print(f"  Сохранено: {OUT}/chr6_windows_Fw_Sw.tsv")
print(f"  Окон с eQTL: {fw_sw['has_eqtl'].sum()} / {len(fw_sw)}")

# Сохраняем per-tissue матрицу
sw.to_csv(f"{OUT}/chr6_Sw_per_tissue.tsv", sep="\t", index=False)
lead.to_csv(f"{OUT}/chr6_lead_eqtl.tsv", sep="\t", index=False)

# Seg_assignments векторно через merge
# Для каждого окна с Fw>0 находим сегмент NIS с максимальным перекрытием

cols_nis = ["chr_win","win_start","win_end",
            "chr_nis","seg_start","seg_end","hap_id","overlap_bp"]
ix = pd.read_csv(f"{OUT}/windows_x_nis.bed", sep="\t",
                 header=None, names=cols_nis)

# Присваиваем win_id
windows = pd.read_csv(f"{OUT}/chr6_windows_sorted.bed", sep="\t",
                      header=None, names=["chr","win_start","win_end"])
windows["win_id"] = range(len(windows))
windows["chr"] = windows["chr"].astype(str)
ix["chr_win"] = ix["chr_win"].astype(str)

ix = ix.merge(
    windows[["chr","win_start","win_end","win_id"]].rename(columns={"chr":"chr_win"}),
    on=["chr_win","win_start","win_end"],
    how="left"
)

# Присваиваем seg_id = индекс строки в NIS (уникальный ID)
# Создаем seg_id из комбинации hap_id + seg_start + seg_end
ix["seg_key"] = ix["hap_id"] + "_" + ix["seg_start"].astype(str) + "_" + ix["seg_end"].astype(str)
seg_map = {k: i for i, k in enumerate(ix["seg_key"].unique())}
ix["seg_id"] = ix["seg_key"].map(seg_map)

# Для каждого окна с Fw>0 берём сегмент с максимальным перекрытием
fw_pos = fw_sw[fw_sw["Fw"] > 0][["win_id"]].copy()
ix_pos = ix[ix["win_id"].isin(fw_pos["win_id"])].copy()
# Для каждого win_id — строка с max overlap_bp
best_seg = (
    ix_pos.loc[ix_pos.groupby("win_id")["overlap_bp"].idxmax()]
    [["win_id","seg_id"]]
    .reset_index(drop=True)
)
best_seg.to_csv(f"{OUT}/chr6_window_seg_ids.tsv", sep="\t", index=False)
print(f"  Сохранены ID сегментов: {OUT}/chr6_window_seg_ids.tsv")
PYEOF

echo "[$(date)] Шаг 5 завершен"

# Аннотация TSS из GENCODE v19
echo "[$(date)] Шаг 6: Извлечение TSS из GENCODE v19"

python3 << 'PYEOF'
import pandas as pd
import numpy as np
import os

WORKDIR = os.path.expanduser("~/nd_pipeline")
OUT = f"{WORKDIR}/results/pipeline_A"
GTF = f"{WORKDIR}/data/gencode/gencode.v19.annotation.gtf"

# Парсим GTF, только строки с feature == "gene" на chr
tss_list = []
with open(GTF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        chrom, source, feature, start, end, score, strand, frame, attrs = parts
        if chrom != "chr6" and chrom != "6":
            continue
        if feature != "gene":
            continue
        start = int(start)
        end = int(end)
        # TSS: для + strand = start, для - strand = end
        tss = start if strand == "+" else end
        # Извлекаем gene_id и gene_name
        gene_id = ""
        gene_name = ""
        for attr in attrs.split(";"):
            attr = attr.strip()
            if attr.startswith("gene_id"):
                gene_id = attr.split('"')[1]
            elif attr.startswith("gene_name"):
                gene_name = attr.split('"')[1]
        tss_list.append({"chr": 6, "tss": tss, "strand": strand,
                         "gene_id": gene_id, "gene_name": gene_name})

tss_df = pd.DataFrame(tss_list)
print(f"  Генов на chr6: {len(tss_df)}")
tss_df.to_csv(f"{OUT}/chr6_tss.tsv", sep="\t", index=False)

# BED для TSS (+-1 bp)
tss_bed = tss_df[["chr","tss","gene_id","gene_name"]].copy()
tss_bed["tss_start"] = tss_bed["tss"] - 1
tss_bed["tss_end"] = tss_bed["tss"]
tss_bed[["chr","tss_start","tss_end","gene_id","gene_name"]].to_csv(
    f"{OUT}/chr6_tss.bed", sep="\t", index=False, header=False)
print(f"  TSS BED сохранён: {OUT}/chr6_tss.bed")
PYEOF
echo "[$(date)] Шаг 6 завершен"

# Расстояние до ближайшего TSS (D_TSS)
echo "[$(date)] Шаг 7: Вычисление D_TSS для каждого окна"
sort -k1,1 -k2,2n "$OUT/chr6_tss.bed" > "$OUT/chr6_tss_sorted.bed"
# chr6_windows.bed (chr, start, end, win_id)
sort -k1,1 -k2,2n "$OUT/chr6_windows.bed" > "$OUT/chr6_windows_with_id_sorted.bed"

# bedtools closest для каждого окна
# -a имеет 4 колонки (chr,start,end,win_id), -b имеет 5 (chr,start,end,gene_id,gene_name)
bedtools closest \
    -a "$OUT/chr6_windows_with_id_sorted.bed" \
    -b "$OUT/chr6_tss_sorted.bed" \
    -d \
    > "$OUT/chr6_windows_dtss.bed"
echo "  D_TSS вычислен: $(wc -l < "$OUT/chr6_windows_dtss.bed") окон"

python3 << 'PYEOF'
import pandas as pd
import numpy as np
import os

WORKDIR = os.path.expanduser("~/nd_pipeline")
OUT = f"{WORKDIR}/results/pipeline_A"
# Читаем рез bedtools closest
cols = ["chr_win","win_start","win_end","win_id",
        "chr_tss","tss_start","tss_end","gene_id","gene_name","D_TSS"]
dtss = pd.read_csv(f"{OUT}/chr6_windows_dtss.bed", sep="\t",
                   header=None, names=cols)

# Категории по D_TSS
def dtss_category(d):
    if d < 5000:
        return "Promoter"
    elif d < 50000:
        return "Near"
    else:
        return "Distal"

dtss["dtss_cat"] = dtss["D_TSS"].apply(dtss_category)
dtss["nearest_gene"] = dtss["gene_name"]

# bedtools closest может вернуть несколько TSS на одинаковом
# оставляем одну запись на win_id (с минимальным D_TSS)
print(f"  Строк до дедупликации: {len(dtss)}")
dtss = dtss.sort_values("D_TSS").drop_duplicates(subset="win_id", keep="first")
print(f"  Строк после дедупликации: {len(dtss)}")

# Объед с Fw/Sw
fw_sw = pd.read_csv(f"{OUT}/chr6_windows_Fw_Sw.tsv", sep="\t")
result = fw_sw.merge(
    dtss[["win_id","D_TSS","dtss_cat","nearest_gene","gene_id"]],
    on="win_id", how="left"
)
result.to_csv(f"{OUT}/chr6_windows_full.tsv", sep="\t", index=False)
print(f"  Итоговая матрица: {OUT}/chr6_windows_full.tsv")
print(f"  Размер: {result.shape}")
print(f"  Категории D_TSS:\n{result['dtss_cat'].value_counts()}")
PYEOF
echo "[$(date)] Шаг 7 завершен"

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
# Для каждого окна [win_start, win_end) HapMap-точки внутри
# среднее rate, если точек нет, интерполируем из ближайших

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
        # Есть HapMap-точки внутри окна — среднее
        recomb_rates[i] = rate_arr[l:r].mean()
    else:
        if l == 0:
            recomb_rates[i] = rate_arr[0]
        elif l >= len(pos_arr):
            recomb_rates[i] = rate_arr[-1]
        else:
            recomb_rates[i] = (rate_arr[l - 1] + rate_arr[l]) / 2.0

# Замена нулей и отриц значений на мин положительное
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
