## ----------------------------READING PACKAGES---------------------------------
# setting random generator seed
set.seed(16748991)
seed_before <- .Random.seed

# creating vector of necessary packages
packages <- c(
  "dplyr",
  "showtext",
  "stringr",
  "data.table",
  "DBI",
  "duckdb",
  "tidyr"
)


## load development version of phiper
library(phiper)

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())

if (any(installed_packages == FALSE)) {
  pak::pkg_install(packages[!installed_packages])
}

# packages loading
invisible(lapply(packages, library, character.only = TRUE))

# add font + phiper use font in all plots
font_add_google("Montserrat", "monte")
phip_use_montserrat()
showtext_auto()

# removing unnecessary variables
rm(list = c('installed_packages', 'packages'))

## ---------------------LOADING DATA--------------------------------------------
# ------------------------------------------------------------------------------
# 1) read data
# ------------------------------------------------------------------------------
data      <- fread("Data/exist.csv")
data_fold <- fread("Data/fold.csv")
# ------------------------------------------------------------------------------
# 2) pivot "exist" to long
# ------------------------------------------------------------------------------
start_idx_exist <- which(startsWith(names(data), "R"))[1]

long_exist <- data |>
  dplyr::select(
    peptide_name,
    dplyr::all_of(names(data)[start_idx_exist:ncol(data)])
  ) |>
  tidyr::pivot_longer(
    cols      = -peptide_name,
    names_to  = "sample_id",
    values_to = "exist"
  ) |>
  dplyr::rename(peptide_id = peptide_name)

# ------------------------------------------------------------------------------
# 3) pivot "fold_change" to long
# ------------------------------------------------------------------------------
start_idx_fc <- which(startsWith(names(data_fold), "R"))[1]

long_fc <- data_fold |>
  dplyr::select(
    peptide_name,
    dplyr::all_of(names(data_fold)[start_idx_fc:ncol(data_fold)])
  ) |>
  tidyr::pivot_longer(
    cols      = -peptide_name,
    names_to  = "sample_id",
    values_to = "fold_change"
  ) |>
  dplyr::rename(peptide_id = peptide_name)

# ------------------------------------------------------------------------------
# 4) join exist + fold_change (by sample_id + peptide_id)
# ------------------------------------------------------------------------------
long_df <- long_exist |>
  dplyr::left_join(
    long_fc,
    by = c("sample_id", "peptide_id")
  )

# ------------------------------------------------------------------------------
# 5) read metadata + clean names + join
# ------------------------------------------------------------------------------
metadata <- fread("Metadata/SNMG_metadata.csv")
colnames(metadata)[1] <- "sample_id"
names(metadata) <- make.unique(names(metadata))

long_df <- long_df |>
  dplyr::left_join(metadata, by = "sample_id")

# replace infinites with max
max_fc <- max(long_df$fold_change[is.finite(long_df$fold_change)], na.rm = TRUE)

idx_inf <- is.infinite(long_df$fold_change)
long_df$fold_change[idx_inf] <- max_fc

## ---------------------------------- DATA EXPORT ------------------------------
## safety-check
data_export <- long_df %>%
  distinct(sample_id, peptide_id, .keep_all = TRUE)

con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:", read_only = FALSE)

# 2. expose the data frame to DuckDB -------------------------------------------
duckdb::duckdb_register(con, "df_tbl", data_export)

# -- df is now a virtual table called "df_tbl" inside the DB
# 3. Persist it to Parquet -----------------------------------------------------
dbExecute(
  con,
  "COPY df_tbl TO 'Data/SNMG.parquet' (FORMAT PARQUET);"
)

# 4. Clean up ------------------------------------------------------------------
dbDisconnect(con, shutdown = TRUE)

rm(list = c("long_df", "data_export", "con", "data", "metadata", "start_idx_exist", "start_idx_fc"))
gc()

