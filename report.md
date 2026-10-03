# Week 1 | Data Cleaning and Preliminary Analysis with R

**Palmer Archipelago Penguins**  
Prepared: 3 October 2026  

## Executive Summary

This report examines the public Palmer Penguins dataset, a compact but realistic ecological dataset with numeric measurements, categorical descriptors, and missing observations. The source contains **344 observations and 8 fields** covering three penguin species, three islands, four physical measurements, sex, and observation year.

The workflow uses base R to inspect structure and summaries, quantify missingness, screen numeric measurements for potential outliers, impute missing values while retaining audit flags, preserve original measurement units while adding min-max scaled copies, encode categorical variables, and explore distributions and relationships. All 344 records remain in the cleaned dataset; it contains no missing cells after preprocessing. The IQR screen flags no observations for the four measured traits, so no measurements were removed as outliers.

Initial patterns are pronounced: Gentoo penguins have the highest average body mass (5,076 g) and flipper length (217.2 mm), while Adelie penguins have the shortest bills on average (38.79 mm). Flipper length and body mass have a strong positive pooled correlation ($r = 0.871$). These are descriptive associations, not causal conclusions; pooling species can confound relationships between measurements.

## Dataset and Provenance

The Palmer Penguins data were assembled from research at Palmer Station, Antarctica, and made available by Palmer Station LTER and Kristen Gorman. The simplified dataset is distributed by the `palmerpenguins` project and is intended for data exploration. The raw CSV used here is downloaded directly from the project's public GitHub repository and included in `data/penguins.csv` for reproducibility.

| Attribute | Description |
|---|---|
| Unit of observation | One adult penguin record |
| Rows and columns | 344 rows x 8 columns |
| Categorical fields | `species`, `island`, `sex` |
| Numeric fields | bill length, bill depth, flipper length, body mass, year |
| Measurement units | Millimetres for bill/flipper; grams for body mass |
| Years | 2007-2009 |
| Missingness | Four measurement fields and sex contain missing values |
| Data licence | CC0, as described by the dataset project |

The four numeric morphology measurements represent bill length, bill depth, flipper length, and body mass. Species composition is unbalanced: Adelie 152, Gentoo 124, and Chinstrap 68. This matters when interpreting pooled means and correlations.

**Source and citation.** Horst, A. M., Hill, A. P., & Gorman, K. B. (2020). *palmerpenguins: Palmer Archipelago (Antarctica) penguin data*. R package version 0.1.0. https://doi.org/10.5281/zenodo.3960218. Original ecological study: Gorman, K. B., Williams, T. D., & Fraser, W. R. (2014). Ecological sexual dimorphism and environmental variability within a community of Antarctic penguins. *PLOS ONE, 9*(3), e90081. https://doi.org/10.1371/journal.pone.0090081.

## Reproducible Workflow

The project uses **R 4.3.3** and base R functions only; no add-on packages are needed. From the project root, run:

```powershell
Rscript analysis.R
```

The script reads `data/penguins.csv` and creates the processed datasets, analysis tables, and figures in `outputs/`. The supplied environment used to prepare these materials was Conda with R 4.3.3; a local R installation can also run the script.

### Initial inspection

The script begins with `read.csv()`, then uses `str()` and `summary()` to inspect field types, ranges, quartiles, and missing values. Reading the CSV directly gives categorical fields as character strings; they are converted to factors after missing categories are handled.

```r
penguins_raw <- read.csv("data/penguins.csv", na.strings = c("NA", ""))
str(penguins_raw)
summary(penguins_raw)
```

**Observed structure:** 344 observations and 8 variables. The four physical measurements are numeric (two stored as integer columns); species, island, and sex are character fields in the source CSV; year is integer. The measurements span 32.1-59.6 mm bill length, 13.1-21.5 mm bill depth, 172-231 mm flipper length, and 2,700-6,300 g body mass.

## Missing-Value Audit and Treatment

There are **19 missing cells** in total. Each of the four morphology variables is missing for 2 records (0.58% per field); sex is missing for 11 records (3.20%). Species, island, and year have no missing values. The complete-case dataset has 333 rows because 11 records have at least one missing value.

Numeric morphology fields are imputed with their observed global median. This is a robust, simple choice for a small preliminary exercise and avoids deleting incomplete records. For every imputed measurement, a companion Boolean `*_was_missing` field records whether the value was filled. Missing sex is assigned an explicit `Unknown` factor level rather than being guessed. No missingness is introduced for fields that were complete.

```r
for (field in numeric_traits) {
  missing_flag <- is.na(penguins_clean[[field]])
  median_value <- median(penguins_clean[[field]], na.rm = TRUE)
  penguins_clean[[paste0(field, "_was_missing")]] <- missing_flag
  penguins_clean[[field]][missing_flag] <- median_value
}
penguins_clean$sex[is.na(penguins_clean$sex)] <- "Unknown"
```

**Result:** 344 records retained; 0 missing cells in the cleaned dataset. The imputed values are convenient for downstream demonstrations, not a claim that the penguins had average measurements. Median imputation can reduce variance and weaken relationships; use the missingness indicators and compare against complete-case results before inferential or predictive use. The two records missing all four measured traits deserve particular attention because those imputations are not supported by other morphology values on the same record.

![Missing-value counts by field](outputs/figure_3_missingness.png)

*Figure 1. Missing records by source field. Each morphology field has two missing observations; sex has eleven.*

## Outlier Screening

For each morphology field, the script computes the standard Tukey fences from observed values: $Q_1 - 1.5 \times IQR$ and $Q_3 + 1.5 \times IQR$. This screening is done before imputation and on the pooled dataset. It identifies **zero flagged values** in all four fields:

| Field | Lower fence | Upper fence | Flagged values |
|---|---:|---:|---:|
| Bill length (mm) | 25.31 | 62.41 | 0 |
| Bill depth (mm) | 10.95 | 23.35 | 0 |
| Flipper length (mm) | 155.50 | 247.50 | 0 |
| Body mass (g) | 1,750 | 6,550 | 0 |

No values were capped or removed. An IQR flag is a review signal, not proof of an error; plausible biological extremes should be retained. Pooled fences can also miss unusual values within a species, so future work should review distributions by species and compare against domain limits.

```r
q1 <- quantile(x, 0.25, na.rm = TRUE)
q3 <- quantile(x, 0.75, na.rm = TRUE)
iqr_value <- q3 - q1
lower_fence <- q1 - 1.5 * iqr_value
upper_fence <- q3 + 1.5 * iqr_value
flagged <- x < lower_fence | x > upper_fence
```

## Transformations and Encoding

### Min-max normalization

For each morphology field, a scaled copy is calculated as $(x - min(x))/(max(x) - min(x))$. Values are therefore on a 0-to-1 scale. The original measurements and their meaningful units remain unchanged in the cleaned file; normalized columns have a `_minmax` suffix. Scaling is useful for algorithms sensitive to feature magnitude, but is not required for descriptive statistics or Pearson correlation. For predictive modeling, the scaling bounds should be learned on the training split only to prevent data leakage.

### Categorical encoding

`species`, `island`, and `sex` are converted to factors, with missing sex represented by `Unknown`. The script creates a full one-hot indicator for each observed category level: three species columns, three island columns, and three sex columns, for **9 indicator columns** total. All levels are retained in the exported design matrix. In a linear model with an intercept, one indicator per categorical field should be selected as a reference or an equivalent contrast scheme used to avoid perfect multicollinearity (the dummy-variable trap).

```r
indicators <- sapply(levels(factor_value), function(level) {
  as.integer(factor_value == level)
})
colnames(indicators) <- paste(field, levels(factor_value), sep = "_")
```

### Deliverable files from the script

| File | Contents |
|---|---|
| `outputs/penguins_clean.csv` | Imputed fields, original units, factor values, and missingness flags |
| `outputs/penguins_scaled.csv` | Cleaned fields plus min-max copies of four traits |
| `outputs/model_matrix.csv` | Scaled numeric predictors and full one-hot indicators |
| `outputs/missingness.csv` | Missing counts and percentages by source field |
| `outputs/outlier_screen.csv` | IQR fences, counts, and source row indices |
| `outputs/descriptive_statistics.csv` | Observed-value statistics for morphology fields |
| `outputs/species_counts.csv` | Source record count by species |
| `outputs/species_means.csv` | Raw observed morphology means by species |
| `outputs/correlations.csv` | Pairwise-complete Pearson correlation matrix |

## Exploratory Analysis

### Descriptive statistics

The following statistics use the original observed measurements, before imputation. `n` excludes missing observations.

| Measurement | n | Mean | Median | SD | Min | Q1 | Q3 | Max |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| Bill length (mm) | 342 | 43.92 | 44.45 | 5.46 | 32.1 | 39.23 | 48.5 | 59.6 |
| Bill depth (mm) | 342 | 17.15 | 17.30 | 1.98 | 13.1 | 15.60 | 18.7 | 21.5 |
| Flipper length (mm) | 342 | 200.92 | 197.00 | 14.06 | 172 | 190 | 213 | 231 |
| Body mass (g) | 342 | 4,201.75 | 4,050 | 801.96 | 2,700 | 3,550 | 4,750 | 6,300 |

```r
summary(penguins_raw)
mean(penguins_raw$body_mass_g, na.rm = TRUE)
sd(penguins_raw$body_mass_g, na.rm = TRUE)
```

![Morphology distributions by species](outputs/figure_1_distributions.png)

*Figure 2. Boxplots of each measured trait by species. Between-species differences are substantial, especially for flipper length and body mass.*

### Species-level patterns

| Species | n | Mean bill length (mm) | Mean bill depth (mm) | Mean flipper (mm) | Mean mass (g) |
|---|---:|---:|---:|---:|---:|
| Adelie | 152 | 38.79 | 18.35 | 190.0 | 3,701 |
| Chinstrap | 68 | 48.83 | 18.42 | 195.8 | 3,733 |
| Gentoo | 124 | 47.50 | 14.98 | 217.2 | 5,076 |

Gentoo penguins have the greatest mean flipper length and body mass in this sample. Chinstrap and Gentoo bill lengths are longer on average than Adelie bill lengths, while Chinstrap bill depth is similar to Adelie. These group averages describe this dataset; they do not establish causes or population-level significance.

### Correlation

Pearson correlations use pairwise-complete observed measurements, so each pair uses available values without first imputing them.

|  | Bill length | Bill depth | Flipper length | Body mass |
|---|---:|---:|---:|---:|
| Bill length | 1.000 | -0.235 | 0.656 | 0.595 |
| Bill depth | -0.235 | 1.000 | -0.584 | -0.472 |
| Flipper length | 0.656 | -0.584 | 1.000 | 0.871 |
| Body mass | 0.595 | -0.472 | 0.871 | 1.000 |

The largest non-self correlation is between flipper length and body mass ($r = 0.871$). Bill length also has a positive association with flipper length ($r = 0.656$), while bill depth is negatively associated with body mass ($r = -0.472$). These pooled correlations may largely reflect species differences, so species-stratified plots and correlations are a sensible next check.

![Bill and flipper measurements](outputs/figure_2_bill_flipper.png)

*Figure 3. Bill length versus flipper length, colored by species; the dark line is a pooled linear trend for orientation only.*

## Key R Output

Selected console results from the executed R script:

```text
Rows: 344  Columns: 8
Complete cases: 333 of 344
Rows retained after imputation: 344
Remaining missing cells in cleaned data: 0
One-hot predictor columns: 9

Species counts: Adelie 152; Chinstrap 68; Gentoo 124
Mean body mass by species: Adelie 3701 g; Chinstrap 3733 g; Gentoo 5076 g
Correlation, flipper length and body mass: 0.871
IQR outlier flags: 0 for each of the four morphology traits
```

Full reproducible outputs are exported as CSV files in `outputs/`; the R script prints the structures, summaries, missingness table, outlier screen, species means, correlations, and record-retention checks to the console.

## Conclusions, Limitations, and Next Steps

This first-pass workflow provides an auditable cleaned dataset and finds clear descriptive differences among species. Missingness is limited but concentrated: two records have all four measured traits missing, and sex is unknown for eleven records. No global IQR outliers were found, so retaining all measurements is justified for this exploratory stage. The body-mass/flipper association is strong in the pooled data, but species composition is a likely confounder.

Limitations include the modest sample size, unequal species counts, imputation uncertainty, the use of pooled outlier fences and correlations, and the absence of inferential tests or sampling-design analysis. Recommended next steps are to stratify distributions and correlations by species, compare complete-case and imputed summaries, investigate missingness patterns against island/year/species, and use training-only scaling if building a predictive model. No causal interpretation is warranted from this exploratory analysis.

## Suggested 32-Hour Learning Plan

The assignment estimates 30-35 hours of effort. This schedule is a suggested learning plan, not a claim about hours spent preparing the report.

| Activity | Hours |
|---|---:|
| Dataset selection, source review, and download | 3 |
| R setup, import, and structure/type inspection | 4 |
| Missingness audit and handling rationale | 5 |
| Outlier screening and validation | 4 |
| Scaling and categorical encoding | 4 |
| Summary statistics and visual exploration | 6 |
| Interpretation, reproducibility checks, and documentation | 6 |
| **Total** | **32** |

## Appendix A. Reproduction Checklist

1. Install R 4.3 or later.
2. Open a terminal in the project root directory.
3. Run `Rscript analysis.R`.
4. Confirm `outputs/` contains the exported CSV tables and three PNG figures.
5. Open `outputs/penguins_clean.csv` and verify there are 344 rows and no missing cells.

The R analysis is in `analysis.R`; the exact downloaded source file is in `data/penguins.csv`.
