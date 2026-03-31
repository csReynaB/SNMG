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
data      <- fread("Data/exist.csv")#[,-"R26P04_66_307032_PCaInnsbruck_Ctrls_A_T_C2"]
data_fold <- fread("Data/fold.csv")#[,-"R26P04_66_307032_PCaInnsbruck_Ctrls_A_T_C2"]
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
metadata <- fread("Metadata/IBD-Berlin_metadata.csv") #%>%
# mutate(
#   group_cmb = ifelse(
#     group_liver != "Controls" & !is.na(Etiology),
#     paste(Etiology, group_sar, sep = "_"),
#     NA
#   )
#   group_test = factor(
#     group_test,
#     levels = c("Cirr+Sar", "Cirr+No Sar", "No cirr+Sar", "Controls")
#   ),
#   Etiology = factor(
#     Etiology,
#     levels = c("HCV", "Alcohol", "MASH", "Others")
#   )
#)

colnames(metadata)[1] <- "sample_id"
names(metadata) <- make.unique(names(metadata))

long_df <- long_df |>
  dplyr::left_join(metadata, by = "sample_id")

# replace infinites with max
max_fc <- max(long_df$fold_change[is.finite(long_df$fold_change)], na.rm = TRUE)

idx_inf <- is.infinite(long_df$fold_change)
long_df$fold_change[idx_inf] <- max_fc

## ---------------------------------- DATA EXPORT ------------------------------
## exporting the data to a .parquet file --> for phiper use

# long_df_subset <- long_df %>%
#   dplyr::filter(!is.na(group_cmb), 
#                 !is.na(Etiology),
#                 !grepl("Others", group_cmb, ignore.case = TRUE),
#                 !grepl("Others", Etiology, ignore.case = TRUE))
# 

## safety-check
data_export <- long_df %>%
  distinct(sample_id, peptide_id, .keep_all = TRUE)

# data_export_subset <- long_df_subset %>%
#   distinct(sample_id, peptide_id, .keep_all = TRUE)


# 1. open (or create) a DuckDB database ----------------------------------------
con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:", read_only = FALSE)

# 2. expose the data frame to DuckDB -------------------------------------------
duckdb::duckdb_register(con, "df_tbl", data_export)

# -- df is now a virtual table called "df_tbl" inside the DB
# 3. Persist it to Parquet -----------------------------------------------------
dbExecute(
  con,
  "COPY df_tbl TO 'Data/IBD-Berlin.parquet' (FORMAT PARQUET);"
)

# 4. Clean up ------------------------------------------------------------------
dbDisconnect(con, shutdown = TRUE)


# 
# # 1. open (or create) a DuckDB database ----------------------------------------
# con <- dbConnect(duckdb::duckdb(), dbdir = ":memory:", read_only = FALSE)
# 
# # 2. expose the data frame to DuckDB -------------------------------------------
# duckdb::duckdb_register(con, "df_tbl", data_export_subset)
# 
# # -- df is now a virtual table called "df_tbl" inside the DB
# # 3. Persist it to Parquet -----------------------------------------------------
# dbExecute(
#   con,
#   "COPY df_tbl TO 'Data/SAR_MUW_subset.parquet' (FORMAT PARQUET);"
# )
# 
# # 4. Clean up ------------------------------------------------------------------
# dbDisconnect(con, shutdown = TRUE)

rm(list = c("long_df", "data_export", "con", "data", "metadata", "start_idx_exist", "start_idx_fc"))
#rm(list = c("long_df", "data_export", "data_export_subset", "con", "data", "metadata", "start_idx_exist", "start_idx_fc"))
gc()

