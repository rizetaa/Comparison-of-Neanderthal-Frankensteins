# Пайплайн eQTL-влияние неандертальской интрогрессии

## пару слов про DAIseg и получение файлов IBS.YRI.grch37.chrN.em.tsv

## Структура файлов
```
~/nd_pipeline/                        # Рабочая директория пайплайна
├── data/
│   ├── raw/                          # IBS.YRI...
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
## 0. Подготовка

### Создание структуры директории
```bash
mkdir -p ~/nd_pipeline/{data/{raw,hapmap,gtex,vindija,gencode},scripts,results,logs}
```
### Копирование tsv-файла на сервер

## 1. Скачивание данных

```bash
cd ~/nd_pipeline
```
### Выполнить все, что указано в 00_download_data.sh

## Установить права на выполнение

```bash
cd ~/nd_pipeline
chmod +x "имя файла"
```
## Тестовый запуск для одной хромосомы
```bash
cd ~/nd_pipeline
CHR=6 bash run_chromosome_without_callable.sh 2>&1 | tee results/chr6/logs/chr6_test.log
```

## 2. Pipeline A (препроцессинг)
```
Подготовка данных для основного анализа 
1. Нарезаем chr6 на окна 1000 bp
2. Вычисляем глобальную частоту интрогрессии Fw для каждого окна
   - Индикатор I=1 если гаплотип покрывает окно >= 50% (>= 500 bp)
   - Fw=(число гаплотипов с I=1) / H
3. Присваиваем бины частот: Zero / Rare / Low / Intermediate / High / Very_High
4. Фильтруем GTEx v7 eQTL: только chr6, FDR < 0.05, биаллельные SNP
5. Пересекаем eQTL с окнами
6. Proxy-clumping: для каждой тройки (окно, ткань, ген) -> lead SNP с max|Z|
7. Агрегируем Sw,t = max_g |Z_lead|
8. Извлекаем TSS из GENCODE v19, вычисляем D_TSS
9. Категоризируем окна: Promoter (<5 kb) / Near (5–50 kb) / Distal (>=50 kb)
10. Вычисляем среднюю скорость рекомбинации из HapMap II
Итоговый файл: results/pipeline_A/chr6_windows_full.tsv
```

## 3. Pipeline B (субпопуляции и поляризация Vindija)
```
Подготовка данных для cпецанализа
1. Группируем гаплотипы по субпопуляциям
2. Вычисляем субпопуляционные частоты F^(k)_w для каждой популяции
3. Пересекаем eQTL с Vindija 33.19 VCF
4. Hard Filtering Vindija:
   - GQ > 30 (качество генотипирования)
   - DP > 10 (глубина покрытия)
   - Исключает транзиции C/T и G/A 
5. Определяем знак D_w,t = sign(beta_neand):
   - +1 если неандертальский аллель повышает экспрессию
   - −1 если понижает
   - NA если аллель неизвестен
6. Включаем только валидные окна, которые достаточно хорошо покрытые callable-областями. Для этого нужно найти долю callable-позиций в окне:
vw = (callable bp in w) / 1000

Итоговый файл: results/pipeline_B/chr6_windows_final.tsv
```

## 4. Основной анализ
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
- Сравнение Sw: интрогрессия (Fw > 0) и контроль (Fw = 0)
- Отдельно для Promoter / Near / Distal
- 10 000 итераций block bootstrap
- Блок = целый NIS-сегмент (для интрогрессии) или геномный участок медианной длины (для контроля)
2.4 Адаптивная интрогрессия
- Кандидаты: Fw > 95 перцентиль, Sw > 95 перцентиль
- Scatter plot Fw * Sw с генами-мишенями
```

## 5. Визуализация
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

## 6. Запуск и графики с callability для всех хромосом

```bash
nohup ./scripts/preprocess_pipeline_A_genome.sh 8 max_absx > 1.log 2>&1 &
nohup ./scripts/main.analysis.genome.sh violin max_absz >2.log 2>&1 &
```

## Запуск для всех

```bash
sudo apt-get install parallel
which parallel
```

### Запустить параллельную обработку

```bash
cd ~/nd_pipeline
bash run_all_chromosomes_without_callable.sh 2>&1 | tee results/logs/all_chromosomes.log
```

### Теперь можно скачать картинки себе локально
