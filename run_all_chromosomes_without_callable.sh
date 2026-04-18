#!/bin/bash
# Скрипт для параллельного запуска пайплайна для всех хромосом
# запуск bash run_all_chromosomes_v2.sh
# GNU parallel

set -euo pipefail

WORKDIR="$HOME/nd_pipeline"
LOG_DIR="$WORKDIR/results/logs"

echo "Параллельный запуск пайплайна для всех хромосом"
echo "Рабочая директория: $WORKDIR"
echo "Директория логов: $LOG_DIR"

# Создаём директорию логов
mkdir -p "$LOG_DIR"

# Список хромосом
CHROMOSOMES=(1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22)
N_PARALLEL=${#CHROMOSOMES[@]}

echo "Хромосомы для обработки: ${CHROMOSOMES[@]}"
echo "Число параллельных задач: $N_PARALLEL"
echo "=========================================="

# Проверяем наличие GNU parallel
if ! command -v parallel &> /dev/null; then
    echo "ОШИБКА: GNU parallel не установлен!"
    echo "Установите: sudo apt-get install parallel"
    exit 1
fi

# Шаг 0: Предобработка GTEx (один раз)
echo ""
echo "Шаг 0: Предобработка GTEx eQTL по хромосомам"

GTEX_BY_CHR="$WORKDIR/data/gtex/by_chr"
if [ -f "$GTEX_BY_CHR/chr1_eqtl_raw.tsv" ] && [ -f "$GTEX_BY_CHR/chr22_eqtl_raw.tsv" ]; then
    echo "  Файлы уже существуют, пропускаем."
    echo "  (Удалите $GTEX_BY_CHR для пересоздания)"
else
    bash "$WORKDIR/scripts/00_prefilter_gtex_v2.sh"
fi

echo ""
echo "Запуск параллельной обработки хромосом"
# Формируем команды для каждой хромосомы
COMMANDS=()
for chr in "${CHROMOSOMES[@]}"; do
    COMMANDS+=("CHR=${chr} bash ${WORKDIR}/run_chromosome_v2.sh")
done

# Запускаем параллельно
echo ""
echo "Запуск параллельной обработки..."

parallel \
    --jobs $N_PARALLEL \
    --keep-order \
    --joblog "${LOG_DIR}/parallel_joblog.tsv" \
    --resume \
    ::: "${COMMANDS[@]}"

# Проверка результатов
echo ""
echo "Проверка результатов"

# Считаем успешно завершённые хромосомы
completed=0
failed=0

for chr in "${CHROMOSOMES[@]}"; do
    chr_results="$WORKDIR/results/chr${chr}"
    
    # Проверяем наличие основных файлов результатов
    if [ -f "$chr_results/pipeline_A/chr${chr}_windows_full.tsv" ] && \
       [ -f "$chr_results/pipeline_B/chr${chr}_windows_final.tsv" ] && \
       [ -f "$chr_results/analysis/bootstrap_results.json" ]; then
        echo "chr${chr}: ✓ Завершено"
        ((completed++))
    else
        echo "chr${chr}: ✗ Ошибка или не завершено"
        ((failed++))
    fi
done

echo "Итого:"
echo "  Завершено: $completed/${#CHROMOSOMES[@]}"
echo "  Ошибок: $failed/${#CHROMOSOMES[@]}"

if [ $failed -eq 0 ]; then
    echo "Все хромосомы обработаны успешно!"
    exit 0
else
    echo "Некоторые хромосомы завершились с ошибками. Проверьте логи."
    exit 1
fi
