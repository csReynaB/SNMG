#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(fs)
  library(quarto)
})

args <- commandArgs(trailingOnly = TRUE)

base_dir <- NULL
group_cols <- NULL
template <- "src/template/phiper_summary_report.qmd"

for (arg in args) {
  if (grepl("=", arg, fixed = TRUE)) {
    parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
    key <- parts[1]
    value <- parts[2]

    if (key == "BASE_DIR") {
      base_dir <- value
    } else if (key == "GROUP_COLS") {
      group_cols <- strsplit(value, ",", fixed = TRUE)[[1]]
      group_cols <- trimws(group_cols)
      group_cols <- group_cols[nzchar(group_cols)]
    } else if (key == "TEMPLATE") {
      template <- value
    }
  }
}

if (is.null(base_dir) || !nzchar(base_dir)) {
  stop("Please provide BASE_DIR=...")
}

if (is.null(group_cols) || length(group_cols) == 0) {
  stop("Please provide GROUP_COLS=group1,group2,...")
}

base_dir <- fs::path_abs(base_dir)
template <- fs::path_abs(template)
project_root <- getwd()
logo_meduni = fs::path(project_root, "assets", "logos", "meduni.svg")
logo_antibody = fs::path(project_root, "assets", "logos", "antibodies.png")

if (!fs::file_exists(template)) {
  stop("Template not found: ", template)
}

message("Base dir: ", base_dir)
message("Template: ", template)
message("Group columns: ", paste(group_cols, collapse = ", "))

for (gc in group_cols) {
  out_dir <- fs::path(base_dir, gc)

  if (!fs::dir_exists(out_dir)) {
    message("Skipping missing group column folder: ", out_dir)
    next
  }

  qmd_copy <- fs::path(out_dir, "phiper_summary_report.qmd")
  fs::file_copy(template, qmd_copy, overwrite = TRUE)

  old_wd <- getwd()
  setwd(out_dir)

  tryCatch({
    quarto::quarto_render(
      input = qmd_copy,
      output_file = paste0("summary_report_", gc, ".html"),
      execute_dir = out_dir,
      execute_params = list(
        base_dir = base_dir,
        group_col = gc,
        top_n = 500,
        tables_open = FALSE,
        include_beta_tables = FALSE,
        include_alpha_tables = FALSE,
        include_tsne3d = TRUE,
        include_pop_tables = TRUE,
        include_pop_plots = TRUE,
        include_delta_tables = TRUE,
        include_delta_feature_plots = TRUE,
        pop_interactive_mode = "link",
        delta_interactive_mode = "embed",
        max_delta_feature_plots = 5,
        logo_meduni = logo_meduni,
        logo_antibody = logo_antibody
      )
    )
    message("Rendered report for: ", gc)
  }, error = function(e) {
    message("Failed for ", gc, ": ", conditionMessage(e))
  }, finally = {
    setwd(old_wd)
    if (fs::file_exists(qmd_copy)) fs::file_delete(qmd_copy)
  })
}
