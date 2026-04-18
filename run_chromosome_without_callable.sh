#!/bin/bash
# Скрипт для запуска полного пайплайна для одной хромосомы
# CHR=N bash run_chromosome_without_callable.sh

set -euo pipefail

CHR=${CHR:-6}
WORKDIR="$HOME/nd_pipeline"
SCRIPTS_DIR="$WORKDIR/scripts"
CHR_RESULTS_DIR="$WORKDIR/results/chr${CHR}"

echo "Запуск пайплайна для chr${CHR}"
echo "Рабочая директория: $WORKDIR"
echo "Директория результатов: $CHR_RESULTS_DIR"

mkdir -p "$CHR_RESULTS_DIR/logs"

echo "Pipeline A: Препроцессинг"
CHR=${CHR} bash "${SCRIPTS_DIR}/01_pipeline_A_preprocess_universal_v2.sh" 2>&1 | tee "${CHR_RESULTS_DIR}/logs/pipeline_A_chr${CHR}.log"

# Проверка успешности
if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo "Ошибка: Pipeline A завершился с ошибкой"
    exit 1
fi

echo "Pipeline A завершен успешно"

echo "Pipeline B: Субпопуляции + поляризация Vindija"
CHR=${CHR} bash "${SCRIPTS_DIR}/02_pipeline_B_subpop_vindija_universal_v2.sh" 2>&1 | tee "${CHR_RESULTS_DIR}/logs/pipeline_B_chr${CHR}.log"

# Проверка успешности
if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo "Ошибка: Pipeline B завершился с ошибкой"
    exit 1
fi

echo "Pipeline B завершен успешно"

echo "Основной анализ"
CHR=${CHR} python3 "${SCRIPTS_DIR}/03_main_analysis_universal_v2.py" 2>&1 | tee "${CHR_RESULTS_DIR}/logs/pipeline_03_chr${CHR}.log"

# Проверка успешности
if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo "ОШИБКА: Основной анализ завершился с ошибкой!"
    exit 1
fi

echo "Основной анализ завершен успешно"

echo "Визуализация"
CHR=${CHR} python3 "${SCRIPTS_DIR}/04_visualize_universal_v2.py" 2>&1 | tee "${CHR_RESULTS_DIR}/logs/pipeline_04_chr${CHR}.log"

# Проверка успешности
if [ ${PIPESTATUS[0]} -ne 0 ]; then
    echo "Ошибка: Визуализация завершилась с ошибкой"
    exit 1
fi

echo "Визуализация завершена успешно"
echo "Пайплайн для chr${CHR} завершен успешно"
echo "Результаты сохранены в: $CHR_RESULTS_DIR"
echo "Логи сохранены в: $CHR_RESULTS_DIR/logs"
