# ------------------------------------------------------------------------------
# Group definitions and ordering
# ------------------------------------------------------------------------------

group_definitions <- list(
  group_test = list(
    group_col = "group_test",
    groups = c(
      "SNMG_with_colocalization", "SNMG_without_colocalization", "Controls"
      )
  ),
  group_SNMG = list(
    group_col = "group_SNMG",
    groups = c(
      "SNMG", "Controls"
    )
  )
)


# ------------------------------------------------------------------------------
# Optional manual comparisons
# Leave as NULL to generate all pairwise contrasts from the active groups
# ------------------------------------------------------------------------------

manual_comparisons <- NULL
manual_longitudinal <- NULL

# Example:
# manual_comparisons <- list(
#   c("Erlangen_BL", "Erlangen_W8"),
#   c("Erlangen_BL", "Erlangen_W12"),
#   c("Erlangen_BL", "Erlangen_W24")
# )
# manual_longitudinal <- c(TRUE, TRUE, TRUE)
