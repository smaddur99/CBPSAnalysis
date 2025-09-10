#' @importFrom magrittr %>%
#' @importFrom rlang .data
#' @importFrom stats setNames
#' @importFrom utils globalVariables
NULL

# Suppress R CMD check notes about global variables used in dplyr/tidyverse functions
utils::globalVariables(c(
  # dplyr/tidyverse variables
  ".", ".data", ".imp", "term", "estimate", "std.error", "statistic", "p.value",
  "conf.low", "conf.high", "se", "ci_lower", "ci_upper", "variable", "row_id",
  "weight", "sample_size", "model_name", "n_miss", "pct_miss", "everything",
  "first", "contains",

  # Custom variables created in data processing
  "bootstrap_significant", "glm_p_value", "glm_estimate", "glm_conf_low",
  "glm_conf_high", "bootstrap_estimate", "bootstrap_conf_low", "bootstrap_conf_high",
  "primary_treatment_p_value", "primary_treatment_estimate", "primary_treatment_significant",
  "n_balanced_vars", "n_total_vars", "sig_symbol", "p_value_formatted", "row_id",
  "pretty_model", "glm_conf_int", "bootstrap_conf_int", "balance_ratio",
  "sample_size_col", "is_significant"
))
