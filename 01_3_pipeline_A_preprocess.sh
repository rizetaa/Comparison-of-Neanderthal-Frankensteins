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

# Вычисление Sw (proxy-clumping)
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

# Для каждого окна с Fw>0 берем сегмент с максимальным перекрытием
fw_pos = fw_sw[fw_sw["Fw"] > 0][["win_id"]].copy()
ix_pos = ix[ix["win_id"].isin(fw_pos["win_id"])].copy()
# Для каждого win_id строка с max overlap_bp
best_seg = (
    ix_pos.loc[ix_pos.groupby("win_id")["overlap_bp"].idxmax()]
    [["win_id","seg_id"]]
    .reset_index(drop=True)
)
best_seg.to_csv(f"{OUT}/chr6_window_seg_ids.tsv", sep="\t", index=False)
print(f"  Сохранены ID сегментов: {OUT}/chr6_window_seg_ids.tsv")
PYEOF

echo "[$(date)] Шаг 5 завершен"

# TSS, D_TSS
# TSS из GENCODE v19
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
        # TSS: для +strand = start, для -strand = end
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
print(f"  TSS BED сохранен: {OUT}/chr6_tss.bed")
PYEOF

echo "[$(date)] Шаг 6 завершен"

# Расстояние до ближайшего TSS (D_TSS)
echo "[$(date)] Шаг 7: Вычисление D_TSS для каждого окна"
sort -k1,1 -k2,2n "$OUT/chr6_tss.bed" > "$OUT/chr6_tss_sorted.bed"
# chr6_windows.bed (chr, start, end, win_id)
sort -k1,1 -k2,2n "$OUT/chr6_windows.bed" > "$OUT/chr6_windows_with_id_sorted.bed"

# bedtools closest для каждого окна
# -a 4 колонки (chr,start,end,win_id), -b 5 (chr,start,end,gene_id,gene_name)
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

# Читаем bedtools closest
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

# bedtools closest возвращает несколько TSS на одинаковом
# Оставляем одну запись на win_id с минимальным D_TSS
print(f"  Строк до дедупликации: {len(dtss)}")
dtss = dtss.sort_values("D_TSS").drop_duplicates(subset="win_id", keep="first")
print(f"  Строк после дедупликации: {len(dtss)}")

# Объединяем с Fw/Sw
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
