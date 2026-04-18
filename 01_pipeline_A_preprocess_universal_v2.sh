#!/bin/bash
# Скрипт 01: Pipeline A, препроцессинг для глобального анализа
# Пререквизиты: bedtools, python (pandas, numpy), tabix
# bash 01_pipeline_A_preprocess_universal_v2.sh

set -euo pipefail

CHR=${CHR:-6}
WORKDIR="$HOME/nd_pipeline"
DATA="$WORKDIR/data"
RES="$WORKDIR/results"
CHR_RES="$RES/chr${CHR}"
NIS_DIR="$WORKDIR/data/raw"
LOG="$CHR_RES/logs/pipeline_A_chr${CHR}.log"

mkdir -p "$CHR_RES/pipeline_A"
mkdir -p "$CHR_RES/logs"
exec > >(tee -a "$LOG") 2>&1

echo "[$(date)] Pipeline A: Старт (chr${CHR})"

WIN_SIZE=1000
NIS="${NIS_DIR}/IBS.YRI.grch37.chr${CHR}.em.tsv"
HAPMAP="$DATA/hapmap/genetic_map_GRCh37_chr${CHR}.txt"
GENCODE_GTF="$DATA/gencode/gencode.v19.annotation.gtf"
GTEX_DIR="$DATA/gtex/GTEx_Analysis_v7_eQTL"
OUT="$CHR_RES/pipeline_A"

# Длины хромосом
declare -A CHR_LEN=(
    [1]=249250621 [2]=243199373 [3]=198022430 [4]=191154276 [5]=180915260
    [6]=171115067 [7]=159138663 [8]=146364022 [9]=141213431 [10]=135534747
    [11]=135006516 [12]=133851895 [13]=115169878 [14]=107349540 [15]=102531392
    [16]=90354753 [17]=81195210 [18]=78077248 [19]=59128983 [20]=63025520
    [21]=48129895 [22]=51304566
)
CHR_LEN_VAL=${CHR_LEN[$CHR]}

echo "[$(date)] Шаг 1: Создание окон 1000 bp для chr${CHR}"
echo -e "${CHR}\t${CHR_LEN_VAL}" > "$OUT/chr${CHR}.genome"

bedtools makewindows \
    -g "$OUT/chr${CHR}.genome" \
    -w $WIN_SIZE \
    > "$OUT/chr${CHR}_windows_1kb.bed"
echo "  Создано окон: $(wc -l < "$OUT/chr${CHR}_windows_1kb.bed")"

echo "[$(date)] Шаг 2: Вычисление глобальной частоты интрогрессии Fw"

# 2a: Конвертируем NIS в BED и считаем H 
python3 << 'PYEOF'
import pandas as pd
import numpy as np
import os

CHR = os.environ.get("CHR", "6")
WORKDIR = os.path.expanduser("~/nd_pipeline")
NIS_DIR = f"{WORKDIR}/data/raw"
OUT = f"{WORKDIR}/results/chr{CHR}/pipeline_A"

# Читаем NIS-файл
nis = pd.read_csv(f"{NIS_DIR}/IBS.YRI.grch37.chr{CHR}.em.tsv", sep=r"\s+")
print(f"  NIS записей: {len(nis)}")
print(f"  Колонки: {list(nis.columns)}")
print(f"  Пример:\n{nis.head(2).to_string()}")

nis["sample_id"] = nis["Sample"].str.rsplit("_", n=1).str[0]
nis["hap"] = nis["Sample"].str.rsplit("_", n=1).str[-1]

# Уникальные гаплотипы
all_haplotypes = nis["Sample"].unique()
H = len(all_haplotypes)
print(f"  Всего гаплотипов H = {H}")

# Сохраняем H
with open(f"{OUT}/H_total.txt", "w") as fh:
    fh.write(str(H))

# Конвертируем NIS в BED (0-based half-open)
nis_bed = nis[["Start", "End", "Sample"]].copy()
nis_bed.insert(0, "chr", str(CHR))
nis_bed["bed_start"] = (nis_bed["Start"] - 1).clip(lower=0)
nis_bed["bed_end"] = nis_bed["End"]
nis_bed[["chr", "bed_start", "bed_end", "Sample"]].to_csv(
    f"{OUT}/nis_segments.bed", sep="\t", index=False, header=False)
print(f"  NIS BED сохранен: {OUT}/nis_segments.bed")

# Медианная длина NIS сегмента
median_seg_len = nis["Length"].median()
print(f"  Медианная длина NIS-сегмента: {median_seg_len:.0f} bp")
with open(f"{OUT}/median_seg_len.txt", "w") as fh:
    fh.write(str(int(median_seg_len)))
PYEOF

# 2b: bedtools intersect окна х NIS-сегменты
echo "  bedtools intersect: окна х NIS-сегменты"

sort -k1,1 -k2,2n "$OUT/chr${CHR}_windows_1kb.bed"  > "$OUT/chr${CHR}_windows_sorted.bed"
sort -k1,1 -k2,2n "$OUT/nis_segments.bed"       > "$OUT/nis_segments_sorted.bed"

bedtools intersect \
    -a "$OUT/chr${CHR}_windows_sorted.bed" \
    -b "$OUT/nis_segments_sorted.bed" \
    -wo \
    > "$OUT/windows_x_nis.bed"

echo "  Пересечений окна х NIS: $(wc -l < "$OUT/windows_x_nis.bed")"

# 2c: Вычисляем Fw векторно
python3 << 'PYEOF'
import pandas as pd
import numpy as np
import os

CHR = os.environ.get("CHR", "6")
WORKDIR = os.path.expanduser("~/nd_pipeline")
OUT = f"{WORKDIR}/results/chr{CHR}/pipeline_A"
WIN_SIZE = 1000
MIN_OVERLAP = WIN_SIZE // 2

# Читаем H
with open(f"{OUT}/H_total.txt") as fh:
    H = int(fh.read().strip())

# Читаем окна
windows = pd.read_csv(f"{OUT}/chr{CHR}_windows_sorted.bed", sep="\t",
                      header=None, names=["chr","win_start","win_end"])
windows["win_id"] = range(len(windows))
windows["chr"] = windows["chr"].astype(str)

# Читаем пересечения
cols = ["chr_win","win_start","win_end",
        "chr_nis","seg_start","seg_end","hap_id","overlap_bp"]
ix = pd.read_csv(f"{OUT}/windows_x_nis.bed", sep="\t",
                 header=None, names=cols)
ix["chr_win"] = ix["chr_win"].astype(str)
print(f"  Пересечений загружено: {len(ix)}")

# Фильтр: перекрытие >= 500 bp
ix_filt = ix[ix["overlap_bp"] >= MIN_OVERLAP].copy()
print(f"  После фильтра (>={MIN_OVERLAP} bp): {len(ix_filt)}")

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

# Объединяем с полным списком окон
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
if "chr" not in fw_df.columns:
    fw_df.insert(0, "chr", int(CHR))
else:
    fw_df["chr"] = int(CHR)

fw_df.to_csv(f"{OUT}/chr{CHR}_windows_Fw.tsv", sep="\t", index=False)
print(f"  Сохранено: {OUT}/chr{CHR}_windows_Fw.tsv")
print(f"  Распределение по бинам:\n{fw_df['freq_bin'].value_counts()}")

# BED для bedtools intersect
fw_df[["chr","win_start","win_end","win_id"]].to_csv(
    f"{OUT}/chr{CHR}_windows.bed", sep="\t", index=False, header=False)
PYEOF

echo "[$(date)] Шаг 2 завершен"

echo "[$(date)] Шаг 3: Обработка GTEx eQTL из предобработанного файла"

EQTL_RAW="$DATA/gtex/by_chr/chr${CHR}_eqtl_raw.tsv"

if [ ! -f "$EQTL_RAW" ]; then
    echo "  Ошибка: Файл $EQTL_RAW не найден"
    exit 1
fi

echo "  Предобработанный файл: $EQTL_RAW ($(du -h "$EQTL_RAW" | cut -f1))"
echo "  Строк: $(wc -l < "$EQTL_RAW")"

python3 << 'PYEOF'
import pandas as pd
import numpy as np
import os
import sys

CHR = os.environ.get("CHR", "6")
WORKDIR = os.path.expanduser("~/nd_pipeline")
OUT = f"{WORKDIR}/results/chr{CHR}/pipeline_A"

raw_path = f"{WORKDIR}/data/gtex/by_chr/chr{CHR}_eqtl_raw.tsv"
tsv_path = f"{OUT}/chr{CHR}_eqtl_all.tsv"
bed_path = f"{OUT}/chr{CHR}_eqtl.bed"

# Читаем уже отфильтрованный файл
eqtl = pd.read_csv(raw_path, sep="\t")
print(f"  Загружено eQTL chr{CHR}: {len(eqtl)}")

if len(eqtl) == 0:
    print(f"  Предупреждение: eQTL chr{CHR} не найдены, создаем пустой BED")
    open(bed_path, "w").close()
    open(tsv_path, "w").close()
    sys.exit(0)
print(f"  Примеры variant_id:\n{eqtl['variant_id'].head(5).to_string()}")

# Парсим variant_id
split_df = eqtl["variant_id"].str.split("_", n=4, expand=True)
split_df.columns = ["vchr","vpos","vref","valt","vbuild"]
eqtl = pd.concat([eqtl, split_df], axis=1)

# Последний фильтр chr
eqtl = eqtl[eqtl["vchr"].isin([f"chr{CHR}", str(CHR)])].copy()
print(f"  eQTL на chr{CHR}: {len(eqtl)}")

# Фильтр на биаллельные SNP
eqtl = eqtl[
    (eqtl["vref"].str.len() == 1) &
    (eqtl["valt"].str.len() == 1)
].copy()
print(f"  После фильтра биаллельных: {len(eqtl)}")
eqtl["vpos"] = eqtl["vpos"].astype(int)

# Z-score
eqtl["Z"] = eqtl["slope"] / eqtl["slope_se"]
eqtl["absZ"] = eqtl["Z"].abs()

eqtl.to_csv(tsv_path, sep="\t", index=False)
print(f"  Сохранен: {tsv_path}")

# BED для пересечения с окнами
eqtl_bed = eqtl[["vchr","vpos","gene_id","tissue","slope","slope_se","Z","absZ","variant_id"]].copy()
eqtl_bed["vstart"] = eqtl_bed["vpos"] - 1
eqtl_bed["vend"]   = eqtl_bed["vpos"]
eqtl_bed["vchr_num"] = str(CHR)
eqtl_bed = eqtl_bed[["vchr_num","vstart","vend","variant_id","gene_id","tissue","slope","slope_se","Z","absZ"]]
eqtl_bed.to_csv(bed_path, sep="\t", index=False, header=False)
print(f"  BED файл eQTL: {bed_path} ({len(eqtl_bed)} записей)")
os.remove(raw_path)
PYEOF

echo "[$(date)] Шаг 3 завершен"

echo "[$(date)] Шаг 4: Пересечение eQTL с окнами"
# Сортируем файлы
sort -k1,1 -k2,2n "$OUT/chr${CHR}_windows.bed" > "$OUT/chr${CHR}_windows_sorted2.bed"
sort -k1,1 -k2,2n "$OUT/chr${CHR}_eqtl.bed"   > "$OUT/chr${CHR}_eqtl_sorted.bed"

# Для каждого eQTL SNP находим окно
bedtools intersect \
    -a "$OUT/chr${CHR}_eqtl_sorted.bed" \
    -b "$OUT/chr${CHR}_windows_sorted2.bed" \
    -wa -wb \
    > "$OUT/chr${CHR}_eqtl_x_windows.bed"

echo "  Пересечений: $(wc -l < "$OUT/chr${CHR}_eqtl_x_windows.bed")"

echo "[$(date)] Шаг 5: Proxy-clumping и вычисление Sw,t"

python3 << 'PYEOF'
import pandas as pd
import numpy as np
import os

CHR = os.environ.get("CHR", "6")
WORKDIR = os.path.expanduser("~/nd_pipeline")
OUT = f"{WORKDIR}/results/chr{CHR}/pipeline_A"

# Результат пересечения
cols = ["chr_eqtl","start_eqtl","end_eqtl","variant_id","gene_id","tissue",
        "slope","slope_se","Z","absZ",
        "chr_win","win_start","win_end","win_id"]

df = pd.read_csv(f"{OUT}/chr{CHR}_eqtl_x_windows.bed", sep="\t",
                 header=None, names=cols)
print(f"  Записей после пересечения: {len(df)}")

# Proxy-clumping для каждой тройки (win_id, tissue, gene_id)
# Выбираем SNP с максимальным |Z|
df["absZ"] = df["Z"].abs()
lead = df.loc[df.groupby(["win_id","tissue","gene_id"])["absZ"].idxmax()].copy()
print(f"  Lead SNP после clumping: {len(lead)}")

# Агрегация: Sw,t = max_g |Z_lead|
sw = lead.groupby(["win_id","tissue"])["absZ"].max().reset_index()
sw.columns = ["win_id","tissue","Sw"]

# Sw = max по всем тканям для каждого окна
sw_max = sw.groupby("win_id")["Sw"].max().reset_index()
sw_max.columns = ["win_id","Sw_max"]

# Объединяем с Fw
fw_df = pd.read_csv(f"{OUT}/chr{CHR}_windows_Fw.tsv", sep="\t")
fw_sw = fw_df.merge(sw_max, on="win_id", how="left")
fw_sw["Sw_max"] = fw_sw["Sw_max"].fillna(0)
fw_sw["has_eqtl"] = (fw_sw["Sw_max"] > 0).astype(int)

fw_sw.to_csv(f"{OUT}/chr{CHR}_windows_Fw_Sw.tsv", sep="\t", index=False)
print(f"  Сохранено: {OUT}/chr{CHR}_windows_Fw_Sw.tsv")
print(f"  Окон с eQTL: {fw_sw['has_eqtl'].sum()} / {len(fw_sw)}")

# Сохраняем per-tissue матрицу
sw.to_csv(f"{OUT}/chr{CHR}_Sw_per_tissue.tsv", sep="\t", index=False)
lead.to_csv(f"{OUT}/chr{CHR}_lead_eqtl.tsv", sep="\t", index=False)

# Seg_assignments, векторно через merge вместо iterrows
cols_nis = ["chr_win","win_start","win_end",
            "chr_nis","seg_start","seg_end","hap_id","overlap_bp"]
ix = pd.read_csv(f"{OUT}/windows_x_nis.bed", sep="\t",
                 header=None, names=cols_nis)

# Присваиваем win_id
windows = pd.read_csv(f"{OUT}/chr{CHR}_windows_sorted.bed", sep="\t",
                      header=None, names=["chr","win_start","win_end"])
windows["win_id"] = range(len(windows))
windows["chr"] = windows["chr"].astype(str)
ix["chr_win"] = ix["chr_win"].astype(str)

ix = ix.merge(
    windows[["chr","win_start","win_end","win_id"]].rename(columns={"chr":"chr_win"}),
    on=["chr_win","win_start","win_end"],
    how="left"
)
# Присваиваем seg_id
ix["seg_key"] = ix["hap_id"] + "_" + ix["seg_start"].astype(str) + "_" + ix["seg_end"].astype(str)
seg_map = {k: i for i, k in enumerate(ix["seg_key"].unique())}
ix["seg_id"] = ix["seg_key"].map(seg_map)

# Для каждого окна с Fw>0 сегмент с максимальным перекрытием
fw_pos = fw_sw[fw_sw["Fw"] > 0][["win_id"]].copy()
ix_pos = ix[ix["win_id"].isin(fw_pos["win_id"])].copy()

best_seg = (
    ix_pos.loc[ix_pos.groupby("win_id")["overlap_bp"].idxmax()]
    [["win_id","seg_id"]]
    .reset_index(drop=True)
)
best_seg.to_csv(f"{OUT}/chr{CHR}_window_seg_ids.tsv", sep="\t", index=False)
print(f"  Сохранены ID сегментов: {OUT}/chr{CHR}_window_seg_ids.tsv")
PYEOF

echo "[$(date)] Шаг 5 завершен"

echo "[$(date)] Шаг 6: Извлечение TSS из GENCODE v19"

python3 << 'PYEOF'
import pandas as pd
import numpy as np
import os

CHR = os.environ.get("CHR", "6")
WORKDIR = os.path.expanduser("~/nd_pipeline")
OUT = f"{WORKDIR}/results/chr{CHR}/pipeline_A"
GTF = f"{WORKDIR}/data/gencode/gencode.v19.annotation.gtf"

# Парсим GTF (только строки с feature == "gene" на chr)
tss_list = []
with open(GTF) as fh:
    for line in fh:
        if line.startswith("#"):
            continue
        parts = line.strip().split("\t")
        if len(parts) < 9:
            continue
        chrom, source, feature, start, end, score, strand, frame, attrs = parts
        if chrom != f"chr{CHR}" and chrom != str(CHR):
            continue
        if feature != "gene":
            continue
            
        start = int(start)
        end = int(end)
        tss = start if strand == "+" else end
        gene_id = ""
        gene_name = ""
        for attr in attrs.split(";"):
            attr = attr.strip()
            if attr.startswith("gene_id"):
                gene_id = attr.split('"')[1]
            elif attr.startswith("gene_name"):
                gene_name = attr.split('"')[1]
        tss_list.append({"chr": int(CHR), "tss": tss, "strand": strand,
                         "gene_id": gene_id, "gene_name": gene_name})

tss_df = pd.DataFrame(tss_list)
print(f"  Генов на chr{CHR}: {len(tss_df)}")
tss_df.to_csv(f"{OUT}/chr{CHR}_tss.tsv", sep="\t", index=False)

# Создаем BED для TSS (+-1 bp)
tss_bed = tss_df[["chr","tss","gene_id","gene_name"]].copy()
tss_bed["tss_start"] = tss_bed["tss"] - 1
tss_bed["tss_end"] = tss_bed["tss"]
tss_bed[["chr","tss_start","tss_end","gene_id","gene_name"]].to_csv(
    f"{OUT}/chr{CHR}_tss.bed", sep="\t", index=False, header=False)
print(f"  TSS BED сохранен: {OUT}/chr{CHR}_tss.bed")
PYEOF

echo "[$(date)] Шаг 6 завершен"

echo "[$(date)] Шаг 7: Вычисление D_TSS для каждого окна"

sort -k1,1 -k2,2n "$OUT/chr${CHR}_tss.bed" > "$OUT/chr${CHR}_tss_sorted.bed"
sort -k1,1 -k2,2n "$OUT/chr${CHR}_windows.bed" > "$OUT/chr${CHR}_windows_with_id_sorted.bed"

# bedtools closest
bedtools closest \
    -a "$OUT/chr${CHR}_windows_with_id_sorted.bed" \
    -b "$OUT/chr${CHR}_tss_sorted.bed" \
    -d \
    > "$OUT/chr${CHR}_windows_dtss.bed"

echo "  D_TSS вычислен: $(wc -l < "$OUT/chr${CHR}_windows_dtss.bed") окон"

python3 << 'PYEOF'
import pandas as pd
import numpy as np
import os

CHR = os.environ.get("CHR", "6")
WORKDIR = os.path.expanduser("~/nd_pipeline")
OUT = f"{WORKDIR}/results/chr{CHR}/pipeline_A"

# Читаем результат bedtools closest
cols = ["chr_win","win_start","win_end","win_id",
        "chr_tss","tss_start","tss_end","gene_id","gene_name","D_TSS"]
dtss = pd.read_csv(f"{OUT}/chr{CHR}_windows_dtss.bed", sep="\t",
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

print(f"  Строк до дедупликации: {len(dtss)}")
dtss = dtss.sort_values("D_TSS").drop_duplicates(subset="win_id", keep="first")
print(f"  Строк после дедупликации: {len(dtss)}")

# Объединяем с Fw/Sw
fw_sw = pd.read_csv(f"{OUT}/chr{CHR}_windows_Fw_Sw.tsv", sep="\t")
result = fw_sw.merge(
    dtss[["win_id","D_TSS","dtss_cat","nearest_gene","gene_id"]],
    on="win_id", how="left"
)
result.to_csv(f"{OUT}/chr{CHR}_windows_full.tsv", sep="\t", index=False)
print(f"  Итоговая матрица: {OUT}/chr{CHR}_windows_full.tsv")
print(f"  Размер: {result.shape}")
print(f"  Категории D_TSS:\n{result['dtss_cat'].value_counts()}")
PYEOF

echo "[$(date)] Шаг 7 завершен"

echo "[$(date)] Шаг 8: Вычисление средней скорости рекомбинации для окон"

python3 << 'PYEOF'
import pandas as pd
import numpy as np
import os

CHR = os.environ.get("CHR", "6")
WORKDIR = os.path.expanduser("~/nd_pipeline")
OUT = f"{WORKDIR}/results/chr{CHR}/pipeline_A"
HAPMAP = f"{WORKDIR}/data/hapmap/genetic_map_GRCh37_chr{CHR}.txt"

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

# Читаем окна
windows = pd.read_csv(f"{OUT}/chr{CHR}_windows_full.tsv", sep="\t")

# Векторный расчет средней скорости рекомбинации
pos_arr  = hapmap["pos"].values
rate_arr = hapmap["rate"].values

win_starts = windows["win_start"].values
win_ends   = windows["win_end"].values

left_idx  = np.searchsorted(pos_arr, win_starts, side="left")
right_idx = np.searchsorted(pos_arr, win_ends,   side="left")
recomb_rates = np.empty(len(windows), dtype=float)

for i in range(len(windows)):
    l, r = left_idx[i], right_idx[i]
    if l < r:
        recomb_rates[i] = rate_arr[l:r].mean()
    else:
        if l == 0:
            recomb_rates[i] = rate_arr[0]
        elif l >= len(pos_arr):
            recomb_rates[i] = rate_arr[-1]
        else:
            recomb_rates[i] = (rate_arr[l - 1] + rate_arr[l]) / 2.0

# Заменяем нули и отрицательные значения
min_rate = recomb_rates[recomb_rates > 0].min() if (recomb_rates > 0).any() else 1e-4
recomb_rates = np.where(recomb_rates <= 0, min_rate, recomb_rates)

windows["recomb_rate"] = recomb_rates

windows.to_csv(f"{OUT}/chr{CHR}_windows_full.tsv", sep="\t", index=False)
print(f"  Обновлена матрица с recomb_rate: {OUT}/chr{CHR}_windows_full.tsv")
print(f"  Средняя скорость рекомбинации: {windows['recomb_rate'].mean():.4f} cM/Mb")
PYEOF

echo "[$(date)] Шаг 8 завершен"
echo "[$(date)] Pipeline A завершен"
echo "  Итоговый файл: $OUT/chr${CHR}_windows_full.tsv"
echo "  Колонки: chr, win_start, win_end, win_id, Fw, n_hap_introgressed,"
echo "           freq_bin, Sw_max, has_eqtl, D_TSS, dtss_cat, nearest_gene,"
echo "           gene_id, recomb_rate"
