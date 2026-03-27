#!/bin/bash
#SBATCH --job-name=nd_pipelineB
#SBATCH --cpus-per-task=16
#SBATCH --mem=40GB
#SBATCH --time=06:00:00
#SBATCH --output=%HOME/nd_pipeline/logs/slurm_pipelineB_%j.log
#SBATCH --error=%HOME/nd_pipeline/logs/slurm_pipelineB_%j.err

echo "[$(date)] Запуск Pipeline B"
mkdir -p "$HOME/nd_pipeline/logs"

# Проверяем результаты Pipeline A
if [ ! -f "$HOME/nd_pipeline/results/pipeline_A/chr6_windows_full.tsv" ]; then
    echo "Результаты Pipeline A не найдены"
    echo "Сначала sbatch run_01_pipeline_A.sh"
    exit 1
fi
bash "$HOME/nd_pipeline/scripts/02_pipeline_B_subpop_vindija.sh"
echo "[$(date)] Pipeline B завершен"
echo "Результаты: $HOME/nd_pipeline/results/pipeline_B/"
