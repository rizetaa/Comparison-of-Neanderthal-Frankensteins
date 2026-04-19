# Пайплайн eQTL-влияние неандертальской интрогрессии

## 1. Получение NIS-сегментов

### Пайплайн для одной хромосомы
```bash
# Ограничиваем данные, выполняем фильтрацию 1000 геномов VCF-файлов
python daiseg.py restrict_1kG -json all.chr22.json -threads 8

# Callability Mask
# Вычисляем геномные окна, которые доступны для анализа (фильтр маски)
python daiseg.py callability -json all.chr22.json -threads 8

# Предварительная обработка
# Объединяем VCF, фильтруем SNP и создает матрицу наблюдений .tsv
python daiseg.py main.prep -json all.chr22.json -threads 8

# Обучение HMM
# Запускаем скрытую марковскую модель для вывода путей интрогрессии
python daiseg.py run -json all.chr22.json

# Используем EM алгоритм для оценки
python daiseg.py run.with.EM -json all.chr22.json
```

### Для всех хромосом

```bash
#!/bin/bash

CONF="all.chr22.json"
THR=8

echo "--- [1/4] Restricting 1000 Genomes ---"
python daiseg.py restrict_1kG -json $CONF -threads $THR

echo "--- [2/4] Calculating Callability ---"
python daiseg.py callability -json $CONF -threads $THR

echo "--- [3/4] Preprocessing Observations (VCF -> TSV) ---"
python daiseg.py main.prep -json $CONF -threads $THR

echo "--- [4/4] Running HMM  ---"
python daiseg.py run -json $CONF

echo "Пайплайн завершился успешно"
```

## Структура файлов
```
~/nd_pipeline/                        # Рабочая директория пайплайна
├── data/
│   ├── raw/                          # NIS-сегменты
│   ├── hapmap/                       # Карты рекомбинации HapMap II
│   ├── gtex/                         # GTEx v7 cis-eQTL
│   ├── vindija/                      # Vindija 33.19 VCF
│   └── gencode/                      # GENCODE v19 GTF
├── scripts/                          # Скрипты пайплайна
│   ├── 01_pipeline_A_without_callable.sh
│   ├── 02_pipeline_B_without_callable.sh
│   ├── 03_main_analysis_without_callable.py
│   ├── 04_visualize_without_callable.py
|   └── ...                           # C callability
├── results/
│   ├── chr1/
│   │   ├── pipeline_A/
│   │   ├── pipeline_B/
│   │   ├── analysis/
│   │   ├── figures/
│   │   └── logs/
│   ├── chr2/
|   ├── ...
|   ├── analysis/genome/max_absz       # C callability
|   ├── pipeline_A/max_absz            # C callability
│   └── ...
└── logs/
```
## 2. Подготовка

### Создаем структуру директории
```bash
mkdir -p ~/nd_pipeline/{data/{raw,hapmap,gtex,vindija,gencode},scripts,results,logs}
```
### Копируем tsv-файлы

## 3. Скачивание данных

```bash
cd ~/nd_pipeline
```
### Выполняем все, что указано в 00_download_data.sh

## Устанавливаем право на выполнение

```bash
cd ~/nd_pipeline
chmod +x "имя файла"
```
## Можно сделать тестовый запуск для одной хромосомы
```bash
cd ~/nd_pipeline
CHR=6 bash run_chromosome_without_callable.sh 2>&1 | tee results/chr6/logs/chr6_test.log
```

## 4. Pipeline A (препроцессинг)
```
Подготовка данных для основного анализа 
1. Нарезаем chr6 на окна 1000 bp
2. Вычисляем глобальную частоту интрогрессии Fw для каждого окна
    Индикатор I=1 если гаплотип покрывает окно >= 50% (>= 500 bp)
    Fw = (число гаплотипов с I=1) / H
3. Присваиваем бины частот: Zero / Rare / Low / Intermediate / High / Very_High
4. Фильтруем GTEx v7 eQTL: только chr6, FDR < 0.05, биаллельные SNP
5. Пересекаем eQTL с окнами
6. Proxy-clumping: для каждой тройки (окно, ткань, ген) -> lead SNP с max|Z|
7. Агрегируем Sw,t = max|Z_lead|
8. Извлекаем TSS из GENCODE v19, вычисляем D_TSS
9. Категоризируем окна: Promoter (<5 kb) / Near (5–50 kb) / Distal (>=50 kb)
10. Вычисляем среднюю скорость рекомбинации из HapMap II
Итоговый файл: results/pipeline_A/chr6_windows_full.tsv
```

## 5. Pipeline B (субпопуляции и поляризация Vindija)
```
Подготовка данных для cпецанализа
1. Группируем гаплотипы по субпопуляциям
2. Вычисляем субпопуляционные частоты F^(k)_w для каждой популяции
3. Пересекаем eQTL с Vindija 33.19 VCF
4. Hard Filtering Vindija:
    GQ > 30 (качество генотипирования)
    DP > 10 (глубина покрытия)
    Исключает транзиции C/T и G/A 
5. Определяем знак D_w,t = sign(beta_neand):
    +1 если неандертальский аллель повышает экспрессию
    −1 если понижает
    NA если аллель неизвестен
6. Включаем только валидные окна, которые достаточно хорошо покрытые callable-областями. Для этого нужно найти долю callable-позиций в окне:
vw = (callable bp in w) / 1000
Итоговый файл: results/pipeline_B/chr6_windows_final.tsv
```

## 6. Основной анализ
```
2.1 Спектр частот (iSFS)
- Гистограмма распределения окон по бинам частот
- Manhattan plot: позиция х Fw
2.2 Двухэтапная модель (очищающий отбор)
Логистическая регрессия
logit(P(Sw > 0)) ~ Fw + log(1 + D_TSS) + log(1 + Rw) + C(chrom)
Линейная регрессия (только окна с eQTL)
Sw ~ Fw + log(1 + D_TSS) + log(1 + Rw) + log(1 + Nw) + C(chrom)
2.3 Валидация: Block Bootstrap Mann-Whitney
 Сравнение Sw: интрогрессия (Fw > 0) и контроль (Fw = 0)
 Отдельно для Promoter / Near / Distal
 10 000 итераций block bootstrap
 Блок = целый NIS-сегмент (для интрогрессии) или геномный участок медианной длины (для контроля)
2.4 Адаптивная интрогрессия
 Кандидаты: Fw > 95 перцентиль, Sw > 95 перцентиль
 Scatter plot Fw и Sw с генами-мишенями
```

## 7. Визуализация
### Для одной хромосомы
```
fig1_iSFS.png/pdf - спектр частот интрогрессии (barplot, log-шкала)
fig2_manhattan.png/pdf - Manhattan plot Fw (chr6)
fig3_violin_boxplot.png/pdf - Violin + Boxplot: Fw_bin и Sw 
fig4_split_violin.png/pdf - Split Violin: интрогрессия и  контроль по D_TSS
fig5_scatter_adaptive.png/pdf - Scatter Fw и Sw
fig6_summary_panel.png/pdf - все результаты
```
### Для всех хромосом
```
adaptive_introgression_scatter.png - кандидаты адаптивной интрогрессии
introgression_bin_barplot.png - спектр частот интрогрессии
introgression_manhattan.png - Manhattan график частоты интрогрессии 
Sw_by_freq_bin_violin.png - Sw по частотному интервалу интрогрессии 
stratified_control_violin.png - контроль Fw = 0 и Fw > 0 c категориями DTSS
```

## 8. Запускаем и строим графики с callability для всех хромосом

```bash
nohup ./scripts/preprocess_pipeline_A_genome.sh 8 max_absx > 1.log 2>&1 &
nohup ./scripts/main.analysis.genome.sh violin max_absz >2.log 2>&1 &
```

### Запускаем для всех

```bash
sudo apt-get install parallel
which parallel
```

### Запускаем параллельную обработку

```bash
cd ~/nd_pipeline
bash run_all_chromosomes_without_callable.sh 2>&1 | tee results/logs/all_chromosomes.log
```

### Теперь можем скачать картинки себе локально
