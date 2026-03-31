# SNMG

Analysis pipeline for **PhIP-Seq data from the SNMG cohort**, built around development workflows from **`phiper`**.

This repository contains code to:

- reshape peptide-level PhIP-Seq matrices into a long-format dataset,
- export the analysis-ready object to **Parquet**,
- run peptide- and taxon-level exploratory and differential analyses with **phiper**,
- generate static and interactive figures,
- render an automated **Quarto HTML summary report**, and
- launch report generation on an **HPC/Slurm** environment.

---

## Repository layout

```text
SNMG/
├── Data/
│   └── SNMG.parquet
├── scripts/
│   └── get_reports.sh
│   └── run_phiper.sh
│   └── submit_phiper_analysis.sh
├── src/
│   ├── R/
│   │   ├── 01-create_phiper_object.R
│   │   ├── 02-run_phiper_analysis.R
│   │   └── 03-render_phiper_reports.R
│   └── template/
│       └── phiper_summary_report.qmd
└── README.md
```

---

## What this repo does

### 1) Build an analysis-ready Parquet file

`src/R/01-create_phiper_object.R` reads peptide-level input matrices and metadata, converts them to long format, joins metadata, replaces infinite fold changes with the maximum finite value, and writes an analysis-ready Parquet file.

Expected inputs:

- `Data/exist.csv`
- `Data/fold.csv`
- `Metadata/SNMG_metadata.csv`

Main output:

- `Data/SNMG.parquet`

This file is designed as the input for downstream `phiper`-based analyses.

### 2) Run the main phiper analysis workflow

`src/R/02-run_phiper_analysis.R` is the main analysis driver. It:

- loads the Parquet data with `phiper::phip_convert()`,
- computes enrichment summaries,
- computes **alpha diversity**,
- computes **beta diversity** using Bray-Curtis distances, PCoA, PERMANOVA, dispersion, and t-SNE,
- runs the **POP framework** for prevalence comparison,
- runs the **DELTA framework** for differential repertoire testing,
- exports tables, RDS objects, SVG plots, and interactive HTML files.

The script is written to support both:

- an **all-samples overview** for a grouping variable, and
- **pairwise comparisons** defined in a project-specific configuration.

### 3) Render an automated HTML report from generated results

`src/R/03-render_phiper_reports.R` copies the Quarto template `src/template/phiper_summary_report.qmd` into each results folder and renders a self-contained summary report per grouping variable.

The report can include:

- enrichment counts,
- alpha diversity plots,
- beta diversity plots,
- interactive 3D t-SNE,
- POP scatterplots,
- DELTA forest plots,
- selected per-feature DELTA plots,
- optional sortable tables.

### 4) Submit report rendering on Slurm

`scripts/get_reports.sh` is a Slurm wrapper that loads `R`, `Quarto`, and `Pandoc`, then runs `src/R/03-render_phiper_reports.R` for one or more grouping columns.

---

## Typical workflow

A typical run looks like this:

```text
Raw / matrix-like inputs
   ├─ Data/exist.csv
   ├─ Data/fold.csv
   └─ Metadata/SNMG_metadata.csv

        ↓

01-create_phiper_object.R
   → Data/SNMG.parquet

        ↓

02-run_phiper_analysis.R
   → results/<group_col>/all/
   → results/<group_col>/<group1>_vs._<group2>/
   → static SVGs
   → interactive HTML plots
   → POP/DELTA tables

        ↓

03-render_phiper_reports.R
   → summary_report_<group_col>.html
```

---

## Requirements

This project uses **R** plus a set of analysis and reporting packages. Based on the scripts currently in the repository, the workflow uses packages including:

- `phiper`
- `dplyr`
- `tidyr`
- `data.table`
- `duckdb`
- `DBI`
- `showtext`
- `ggplot2`
- `ggpubr`
- `mgcv`
- `openxlsx`
- `plotly`
- `htmlwidgets`
- `svglite`
- `quarto`
- `fs`
- `readr`
- `readxl`
- `knitr`
- `htmltools`

On HPC, the provided Slurm script expects modules similar to:

- `R`
- `Quarto`
- `Pandoc/3.9`

---

## Input data format

### Matrix inputs

The object-building script assumes that `exist.csv` and `fold.csv` contain:

- one peptide identifier column named `peptide_name`, and
- sample columns starting with `R...`

The script pivots those wide matrices to long format and renames:

- `peptide_name` → `peptide_id`
- matrix column names → `sample_id`

### Metadata

The metadata file is expected at:

- `Metadata/SNMG_metadata.csv`

The first column is renamed to `sample_id` before joining to the long-format peptide table.

---

## Running the pipeline

### Step 1 — create the Parquet file

Run from the repository root:

```bash
Rscript src/R/01-create_phiper_object.R
```

This creates:

```bash
Data/SNMG.parquet
```

### Step 2 — run the main phiper analysis

```bash
Rscript src/R/02-run_phiper_analysis.R \
  PROJECT_DIR=SNMG \
  PARQUET_NAME=SNMG.parquet \
  ACTIVE_GROUP=group_test \
  N_CORES=30 \
  MAX_GB=40
```

Useful command-line parameters supported by the script include:

- `PROJECT_DIR`
- `PARQUET_NAME`
- `ACTIVE_GROUP`
- `DEFAULT_LONGITUDINAL`
- `N_CORES`
- `MAX_GB`
- `LOG`
- `LOG_FILE`
- `FORCE`
- `ALL`

---

## Results structure

The analysis script writes outputs under:

```text
<PROJECT_DIR>/results/<group_col>/
```

For each grouping variable, outputs typically include:

### `all/`

Global summaries across all levels of a grouping column, such as:

- `enrichment_counts.svg`
- `alpha_diversity/`
- `beta_diversity/`

### `<group1>_vs._<group2>/`

Pairwise comparison outputs, including:

- `enrichment_counts.svg`
- `alpha_diversity/`
- `beta_diversity/`
- `POP_framework/`
- `DELTA_framework/`
- comparison-specific `.rds`, `.csv`, `.xlsx`, `.svg`, and `.html` files

---

## Automated report generation

The Quarto template at:

```bash
src/template/phiper_summary_report.qmd
```

collects existing result files and builds an HTML report directly inside the selected results folder.

The rendered report summarizes:

- analysis overview,
- detected comparisons,
- alpha diversity,
- beta diversity,
- interactive t-SNE,
- POP results,
- DELTA results,
- selected feature-level plots.

### Render reports manually

```bash
Rscript src/R/03-render_phiper_reports.R \
  BASE_DIR=/absolute/path/to/SNMG/results \
  GROUP_COLS=group_test
```

You can render multiple grouping columns at once with a comma-separated list:

```bash
Rscript src/R/03-render_phiper_reports.R \
  BASE_DIR=/absolute/path/to/SNMG/results \
  GROUP_COLS=group_test,group_SNMG
```

Output files are named like:

```bash
summary_report_<group_col>.html
```

---

## Running report generation on Slurm

Submit the provided helper script from the repository root:

```bash
sbatch scripts/get_reports.sh /absolute/path/to/SNMG/results group_test
```

Or for several grouping variables:

```bash
sbatch scripts/get_reports.sh /absolute/path/to/SNMG/results group_test,group_SNMG
```

The Slurm job currently uses:

- job name: `SNMG-report`
- time: `00:20:00`
- CPUs: `1`
- memory: `10G`

Logs are written to:

```bash
logs/%x_%j.out
logs/%x_%j.err
```

---

## Notes on project-specific configuration

The main analysis script expects a project-specific configuration file at:

```bash
<PROJECT_DIR>/R/group_config.R
```

This configuration is used to define:

- active grouping variables,
- group levels,
- pairwise comparisons,
- whether comparisons are longitudinal,
- palettes and related settings.

If you reuse this repository for another cohort or project, this is one of the main places to adapt.

---

## Outputs included in this repository

This repository already contains:

- analysis scripts using `phiper`,
- a Parquet dataset in `Data/`,
- Quarto report-generation code,
- a Slurm submission helper,
- assets for report branding/logos.

---


## Acknowledgments

Report template and analysis workflow are built around ongoing development of **phiper** and downstream PhIP-Seq reporting utilities.

---

## Contact

Maintained by **Carlos S. Reyna-Blanco**.
