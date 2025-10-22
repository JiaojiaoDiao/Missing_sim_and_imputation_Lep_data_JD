# 🦋  Missing_sim_and_imputation_Lep_data_JD (v2025-10-20)

**Description:**
This revised R script performs trait imputation and uncertainty assessment for Lepidoptera species using the original, non-dummy-encoded trait data. It evaluates three imputation methods (TaxPaint, MissForest with and without phylogeny) under multiple random seeds using Macro-F1 scores to quantify categorical uncertainty.

---

## 🧬 Overview

The pipeline simulates missingness, applies multiple imputation methods, and evaluates performance across five replicate seeds to estimate model stability and uncertainty.
It is fully compatible with the **ALUS_GLOBI Lepidoptera Imputation Framework**.

---

## ⚙️ Workflow Summary

Each replicate (`Seed 1`–`Seed 5`) follows the same four-step structure:

1. **Data Preprocessing**

   * Load the complete Lepidoptera dataset (`Complete_Leptraits_for_analysis_no_dummy.tsv`)
   * Prune the reference phylogenetic tree (`Leptree.newick`)
   * Simulate missingness (MCAR, 5–50%)
   * Convert categorical variables to factors (`convert_all_numeric_to_factor`)
   * Generate predictor lists and append phylogenetic eigenvectors

2. **CASE 1 – TaxPaint**

   * Impute missing traits using genus-level means (`AverageLepTraits_new.txt`)
   * Resolve logical conflicts in trait groups (e.g., *specialist*, *generalist*, *mixed*)
   * Evaluate categorical accuracy using **Macro-F1 scores**

3. **CASE 2 – MissForest (No Phylogeny)**

   * Perform random forest-based imputation on trait data
   * Evaluate Macro-F1 performance using `calculate_macro_f1_missforest()`

4. **CASE 3 – MissForest (With Phylogeny)**

   * Incorporate phylogenetic eigenvectors into MissForest imputations
   * Evaluate with the same Macro-F1 metrics

Each seed run produces a new set of results, ensuring robustness and allowing assessment of stochastic variation in imputation accuracy.

---

## 📂 Output Structure

All outputs are saved automatically to the `Result/` directory:

```
Result/
├─ Seed_1_method1.csv   # TaxPaint results (Macro-F1)
├─ Seed_1_method2.csv   # MissForest (no phylogeny)
├─ Seed_1_method3.csv   # MissForest (with phylogeny)
├─ Seed_2_method1.csv
├─ ...
└─ Seed_5_method3.csv
```

Each file includes Macro-F1 results for multiple missingness levels (5–50%), enabling comparative evaluation across seeds and methods.

---

## 📊 Key Revisions

| Feature              | Description                                        | Benefit                                       |
| -------------------- | -------------------------------------------------- | --------------------------------------------- |
| **Metric update**    | Introduced **Macro-F1** for categorical evaluation | Captures both precision & recall balance      |
| **Data source**      | Uses original (non-dummy) dataset                  | Maintains factor integrity                    |
| **Replication**      | 5 replicate seeds (1–5)                            | Quantifies uncertainty & improves reliability |
| **Method coverage**  | TaxPaint, MissForest (± phylogeny)                 | Unified testing framework                     |
| **Factor handling**  | Automatic conversion of all categorical variables  | Avoids dummy-variable distortion              |
| **Automated export** | Saves CSV results per method & seed                | Facilitates reproducible aggregation          |

---

## 🧠 Methodological Notes

* **Missingness Simulation:**
  Conducted under MCAR (Missing Completely at Random) conditions with missing fractions from 5% to 50%.

* **Evaluation Metric:**
  Macro-averaged F1 score computed across all categorical traits and missingness levels:
  [
  F1 = 2 \times \frac{Precision \times Recall}{Precision + Recall}
  ]

* **Phylogenetic Integration:**
  Phylogenetic eigenvectors derived from a pruned `Leptree.newick` tree using `AppendEigenvectorsPVR_List()`.

* **Uncertainty Quantification:**
  Repetition across five random seeds allows computation of the mean and variance of F1 scores for each method.

---

## 🧩 Dependencies

**R packages used:**
`tidyverse`, `dplyr`, `missForest`, `ape`, `picante`, `PVR`, `TDIP`, `writexl`, `purrr`, and others (see code header for full list).
Custom functions are sourced from:

```
Functions/
├─ DataHandling_Functions.R
├─ Imputation_Functions.R
├─ Phylo_Functions.R
├─ Simpute_Functions.R
└─ Additional_Functions.R
```

---

## 🧾 Citation

If you use this workflow, please cite the corresponding data and method sources:

* **Ries et al. (2022)** – LepTraits dataset
* **Kawahara et al. (2023)** – Feeding and specialization traits
* **TaxPaint** (Matgend, GitHub)
* **MissForest** (Stekhoven & Bühlmann, 2012, *Bioinformatics*)

---

## 📅 Version History

* **v2025-10-20 (JD revision)**

  * Added Macro-F1 scoring
  * Replaced dummy data with original trait dataset
  * Added five replicate seed runs for uncertainty estimation

---

## ✅ Summary

This revision converts the original exploratory workflow into a **fully automated, uncertainty-aware imputation framework** for categorical Lepidoptera traits.
It ensures reproducibility, phylogenetic consistency, and robust evaluation across imputation methods and random seeds.

