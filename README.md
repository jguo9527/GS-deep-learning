# Project: GS-deep-learning

This folder document analysis scripts for genomic selection (GS) experiments and machine-learning models that combine genotype and phenotype data.

---

**Contents**
- `GS-deep-learning/`: R scripts for preprocessing phenotype/genotype data, genomic prediction, and neural-network based models.
- For more detials, visit publication: Multi-Trait Genomic Prediction of Yield-Related Traits in US Soft Wheat under Variable Water Regimes (https://www.mdpi.com/2073-4425/11/11/1270)
---

**Overview**

The scripts in `GS-deep-learning` provide a mostly end-to-end analysis flow:

- `cv_BGLR_v_DL.R` — Main pipeline:
  - Load and combine multi-year phenotype CSVs, calculate adjusted means/BLUPs with `lmerTest`.
  - Import and impute genotype data (HapMap -> MARKOV imputation via `NAM`).
  - Build hybrid genotypes, PCA/clustering (adegenet) and produce PCA plots.
  - Cluster training and validation sets by genetic grouping (Discriminant Analysis)
  - Run genomic prediction with `BGLR` (RKHS), `rrBLUP`, and example `tensorflow`/`keras` ML sections.
  - Produces prediction CSVs and PNG figures.

- `cv_ML_MMDL.R` — Multi-output deep-learning (MMDL): imputation, design-matrix construction, Keras multi-output network, CV loops and prediction exports.
- `cv_ML_MSDL.R` — Single-output / simplified deep-learning (SMDL): single-trait training loops, GPU toggle, predictions and performance summaries.
- `rsm.R` — Response surface modeling (`rsm`) to explore tuning surfaces (epochs × neurons) and generate contour/perspective plots.

---

**Inputs & expected files**
- Phenotype CSVs (author used files named like `g2f_2014_hybrid_data_clean.csv`).
- Genotype input (HapMap or marker table) and optional preprocessed marker tables.
- Many scripts use `file.choose()`; replace with explicit paths for batch runs.

**Typical outputs**
- Prediction tables (CSV): e.g., `predy_GBLUP_*.csv`, `predy_rrblup_*.csv`, `predy_MMDL_*.csv`.
- Diagnostic plots (PNG): PCA, MDS, training curves, and response surface visualizations.

---

**Key R packages**
- `lmerTest`, `tidyr`, `stringr`, `NAM`, `adegenet`, `ggplot2`, `BGLR`, `rrBLUP`, `tensorflow`, `keras`, `mice`, `beepr`, `rsm`.

Install examples:

```bash
R -e "install.packages(c('tidyr','stringr','ggplot2','data.table','gplots','rsm'))"
R -e "if (!requireNamespace('BiocManager', quietly=TRUE)) install.packages('BiocManager'); BiocManager::install(c('adegenet'))"
```

For `tensorflow`/`keras` follow the official installation guide: https://tensorflow.rstudio.com/.

---
