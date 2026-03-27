#!/bin/bash
#SBATCH --job-name=nd_analysis
#SBATCH --cpus-per-task=30
#SBATCH --mem=80GB
#SBATCH --time=08:00:00
#SBATCH --output=%HOME/nd_pipeline/logs/slurm_analysis_%j.log
#SBATCH --error=%HOME/nd_pipeline/logs/slurm_analysis_%j.err

echo "[$(date)] Запуск основного анализа"
mkdir -p "$HOME/nd_pipeline/logs"

# Проверяем результаты Pipeline A
if [ ! -f "$HOME/nd_pipeline/results/pipeline_A/chr6_windows_full.tsv" ]; then
    echo "Результаты Pipeline A не найдены"
    echo "Сначала sbatch run_01_pipeline_A.sh"
    exit 1
fi
python3 "$HOME/nd_pipeline/scripts/03_main_analysis.py"
echo "[$(date)] Основной анализ завершен"
echo "Результаты: $HOME/nd_pipeline/results/analysis/"
