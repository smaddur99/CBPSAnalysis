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
  "n_balanced_vars", "n_total_vars", "sig_symbol", "p_value_formatted",
  "pretty_model", "glm_conf_int", "bootstrap_conf_int", "balance_ratio",
  "sample_size_col", "is_significant",

  # Model weights column used in glm()/polr() calls
  ".w",

  # Pooled results and odds ratio columns
  "fmi", "exp_estimate", "exp_ci_lower", "exp_ci_upper",

  # Balance summaries (npcbps_weighted_analysis, extract_cbps_results)
  "covariate", "covariate_type", "metric", "before", "after", "vr_after",
  "mean_abs_after", "max_abs_after", "min_vr_after", "max_vr_after", "balance_status",
  "type", "balance_stat", "balance_stat_unadj", "var_ratio_adj",
  "max_balance_stat", "min_var_ratio", "max_var_ratio"
))
