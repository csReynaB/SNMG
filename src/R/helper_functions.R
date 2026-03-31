# ------------------------------------------------------------------------------
# Build active configuration
# ------------------------------------------------------------------------------

build_group_config <- function(active_group_name,
                               default_longitudinal = FALSE,
                               group_definitions,
                               fallback_palette = phip_palette,
                               manual_comparisons = NULL,
                               manual_longitudinal = NULL) {
  if (!active_group_name %in% names(group_definitions)) {
    stop(
      "Unknown active_group_name: ", active_group_name,
      ". Available: ", paste(names(group_definitions), collapse = ", ")
    )
  }
  
  group_palette_full <- make_group_palette_full(
    group_definitions = group_definitions,
    fallback_palette = fallback_palette
  )
  
  cfg <- group_definitions[[active_group_name]]
  group_col <- cfg$group_col
  groups <- cfg$groups
  
  if (length(groups) == 0) {
    stop("No groups defined for active_group_name = ", active_group_name)
  }
  
  missing_labels <- setdiff(groups, names(group_palette_full))
  if (length(missing_labels) > 0) {
    stop(
      "These labels are missing from group_palette_full: ",
      paste(missing_labels, collapse = ", ")
    )
  }
  
  if (is.null(manual_comparisons)) {
    comparisons <- combn(groups, 2, simplify = FALSE)
    longitudinal <- rep(default_longitudinal, length(comparisons))
  } else {
    cmp_labels <- unique(unlist(manual_comparisons))
    bad_labels <- setdiff(cmp_labels, groups)
    
    if (length(bad_labels) > 0) {
      stop(
        "These labels in manual_comparisons are not present in active group '",
        active_group_name, "': ",
        paste(bad_labels, collapse = ", ")
      )
    }
    
    comparisons <- manual_comparisons
    
    if (is.null(manual_longitudinal)) {
      longitudinal <- rep(default_longitudinal, length(comparisons))
    } else {
      if (length(manual_longitudinal) != length(comparisons)) {
        stop(
          "manual_longitudinal must have the same length as manual_comparisons."
        )
      }
      longitudinal <- manual_longitudinal
    }
  }
  
  list(
    active_group_name = active_group_name,
    default_longitudinal = default_longitudinal,
    group_col = group_col,
    groups = groups,
    comparisons = comparisons,
    longitudinal = longitudinal,
    group_palette = group_palette_full[groups],
    group_palette_full = group_palette_full
  )
}

# ------------------------------------------------------------------------------
# Build global palette from all labels across all group definitions
# ------------------------------------------------------------------------------
make_group_palette_full <- function(group_definitions,
                                    fallback_palette = phip_palette) {
  all_labels <- unique(unlist(lapply(group_definitions, `[[`, "groups")))
  n_labels <- length(all_labels)
  if (n_labels >= 12) {
    palette_idx <- c(seq_len(n_labels)[-12], n_labels + 1)
  } else {
    palette_idx <- seq_len(n_labels)
  }
  
  if (max(palette_idx) > length(fallback_palette)) {
    stop(
      "Not enough colors in palette: requested color index ",
      max(palette_idx), " but palette has only ",
      length(fallback_palette), " colors. Extend phip_palette."
    )
  }
  
  stats::setNames(
    fallback_palette[seq_along(all_labels)],
    all_labels
  )
}


# ------------------------------------------------------------------------------
# Other helpers for phiper analysis
# ------------------------------------------------------------------------------


# columns that are always retained when exporting data
base_cols <- c("sample_id", "peptide_id", "group_char", "exist")
special_features <- c(
    "is_IEDB_or_cntrl", "is_auto", "is_infect", "is_EBV",
    "is_toxin", "is_PNP", "is_EM", "is_MPA", "is_patho",
    "is_probio", "is_IgA", "is_flagellum", "signalp6_slow",
    "is_topgraph_new", "is_allergens", "anno_is_fungi", "anno_is_food",
    "anno_is_homo_sapiens", "anno_is_lacto_phage"
)
  
# Function to calculate Kulczynski similarity ----------------------------------
kulczynski_mat <- function(A, B) {
  a <- t(A) %*% B   # (nA x nB) shared presences
  nA <- colSums(A)  # length nA
  nB <- colSums(B)  # length nB
  termA <- sweep(a, 1, nA, "/")  # divide each row by nA
  termB <- sweep(a, 2, nB, "/")  # divide each col by nB
  kulc <- 0.5 * (termA + termB)
  
  # If a sample has nA==0 or nB==0, define similarity as 0
  kulc[!is.finite(kulc)] <- 0
  kulc
}

# Function to format p-values nicely -------------------------------------------
format_pval <- function(p, alpha = 0.05) {
  
  # helper to drop trailing zeros (e.g. "1.00" → "1", "0.500" → "0.5")
  drop_zeros <- function(x) sub("\\.?0+$", "", x)
  
  if (is.na(p)) return("NA")
  
  # -------------------------------------------------------------------------
  # 1. Non-significant (p > alpha)
  # -------------------------------------------------------------------------
  if (p > alpha) {
    raw <- formatC(p, digits = 2, format = "f")
    raw <- drop_zeros(raw)
    return(paste0("ns [", raw, "]"))
  }
  
  # -------------------------------------------------------------------------
  # 2. Normal fixed-decimal formatting (0.001 ≤ p ≤ alpha)
  # -------------------------------------------------------------------------
  if (p >= 0.001) {
    raw <- formatC(p, digits = 3, format = "f")
    raw <- drop_zeros(raw)
    return(raw)
  }
  
  # -------------------------------------------------------------------------
  # 3. Scientific notation (< 0.001)
  # -------------------------------------------------------------------------
  raw <- formatC(p, digits = 2, format = "e")
  
  # remove unnecessary zeros: "1.00e-05" → "1e-05"
  raw <- sub("([0-9]+)\\.0+e", "\\1e", raw)
  
  # remove trailing zeros inside "1.10e-04" → "1.1e-04"
  raw <- sub("([0-9]+\\.[0-9]*[1-9])0+e", "\\1e", raw)
  
  return(raw)
}


# keep <- ps_cmp %>%
#   group_by(group_char, peptide_id) %>%
#   summarise(prevalence = mean(exist), .groups = "drop") %>%
#   group_by(peptide_id) %>%
#   filter(any(prevalence >= 0.05)) %>%
#   pull(peptide_id) %>%
#   unique()


# ------------------------------------------------------------------------------
# I/O helpers
# ------------------------------------------------------------------------------
# save an R object to an RDS file, creating parent directories if needed
save_rds_safe <- function(x, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  saveRDS(x, file = path)
}

# select requested columns (base + extras), optionally drop rows with NA in
# selected columns, collect to R, and save to RDS
make_and_save <- function(data, out_path, extra_vars = NULL, drop_na = TRUE) {
  cols_to_select <- unique(c(base_cols, extra_vars %||% character(0)))
  
  avail_cols <- colnames(data$data_long)
  missing_cols <- setdiff(cols_to_select, avail_cols)
  if (length(missing_cols) > 0) {
    message("Skipping missing columns: ", paste(missing_cols, collapse = ", "))
  }
  
  df <- data %>%
    dplyr::select(dplyr::any_of(cols_to_select))
  
  if (isTRUE(drop_na)) {
    cols_present <- intersect(cols_to_select, colnames(df$data_long))
    df <- df %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(cols_present), ~ !is.na(.x)))
  }
  
  df <- df %>% dplyr::collect()
  save_rds_safe(df, out_path)
  invisible(df)
}


  extract_tbl <- function(obj) {
    if (is.data.frame(obj)) {
      return(tibble::as_tibble(obj))
    }
    for (nm in c("data", "table", "tbl", "df", "result", "results")) {
      if (!is.null(obj[[nm]])) {
        return(tibble::as_tibble(obj[[nm]]))
      }
    }
    out <- try(tibble::as_tibble(obj), silent = TRUE)
    if (!inherits(out, "try-error")) {
      return(out)
    }
    out <- try(as.data.frame(obj), silent = TRUE)
    if (!inherits(out, "try-error")) {
      return(tibble::as_tibble(out))
    }
    stop("Cannot extract a data table from the POP result object.")
  }


    
  downsample_for_static <- function(df, prop = 0.1, seed = 1L) {
    if (is.null(df) || !nrow(df)) {
      return(df)
    }
    n <- nrow(df)
    size <- max(1L, floor(n * prop))
    if (size >= n) {
      return(df)
    }
    set.seed(seed)
    dplyr::slice_sample(df, n = size)
  }



  get_binary_and_ids <- function(feature, peplib, tax_cols,
                                 peptide_col = "peptide_id") {
    
    if (feature %in% special_features) {
      vals <- peplib[[feature]]
      present <- !is.na(vals) & as.logical(vals)
      
    } else if (feature %in% names(peplib)) {
      vals <- peplib[[feature]]
      present <- as.logical(vals)
      present[is.na(present)] <- FALSE
      
    } else {
      if (length(tax_cols) == 0L) {
        stop("No taxonomic columns found in peptide library.")
      }
      matches <- lapply(tax_cols, function(col) {
        vals <- peplib[[col]]
        !is.na(vals) & vals == feature
      })
      present <- Reduce(`|`, matches)
    }
    
    peptide_ids <- as.character(peplib[[peptide_col]][present])
    peptide_ids <- unique(peptide_ids)
    list(present = present, peptide_ids = peptide_ids)
  }


  add_background_static <- function(p, bg,
                                    size = 0.8,
                                    alpha = 0.12,
                                    color = "#808080") {
    if (is.null(bg) || !nrow(bg)) return(p)
    
    bg_layer <- ggplot2::geom_point(
      data = bg,
      mapping = ggplot2::aes(x = percent1, y = percent2),
      inherit.aes = FALSE,
      color = color,
      size = size,
      alpha = alpha,
      show.legend = FALSE
    )
    
    p$layers <- c(list(bg_layer), p$layers)
    p
  }