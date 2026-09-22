#' Extract NPCBPS Results into Structured Format
#'
#' Extracts key results from \code{npcbps_weighted_analysis()} into a structured
#' format suitable for tables, visualization, and further analysis.
#'
#' Balance is summarized across ALL imputations using the statistic appropriate
#' to the treatment type:
#' \itemize{
#'   \item Binary treatment: |SMD| (primary; raw difference for binary covariates)
#'     plus variance ratio (secondary)
#'   \item Continuous treatment: |r|, the weighted treatment-covariate correlation
#' }
#' Single-imputation GLM estimates are reported with robust (sandwich) standard
#' errors. Pooled (Rubin's rules) estimates remain the primary results.
#'
#' @param cbps_results List. Output from \code{npcbps_weighted_analysis()}
#' @param include_balance_details Logical. Return full list (TRUE) or a simple
#'   estimates data frame (FALSE). Default TRUE.
#' @param smd_threshold Numeric. |SMD| threshold for binary treatments (default 0.1)
#' @param corr_threshold Numeric. |r| threshold for continuous treatments (default 0.1)
#' @param vr_range Numeric length 2. Acceptable variance ratio range (default c(0.5, 2))
#' @param robust_type Character. Sandwich estimator type (default "HC0")
#'
#' @return If \code{include_balance_details = TRUE}, a list containing:
#' \itemize{
#'   \item \code{estimates}: Imputation-1 GLM estimates (robust SEs) + bootstrap
#'   \item \code{balance}: Per-covariate balance summarized across imputations
#'     (mean and max of the balance statistic, VR range, balance status).
#'     Use this for the supplementary balance tables.
#'   \item \code{balance_by_imputation}: Long table, one row per covariate per imputation
#'   \item \code{sample_sizes}: Sample size information
#'   \item \code{summary_stats}: Overall balance metrics across imputations
#'   \item \code{pooled_estimates}: Rubin's rules pooled estimates (PRIMARY RESULTS)
#'   \item \code{treatment_type}: "binary" or "continuous"
#'   \item \code{balance_statistic}: "|SMD|" or "|r|"
#' }
#'
#' @export
#' @importFrom dplyr mutate rename select group_by summarise ungroup case_when left_join any_of across matches n
#' @importFrom tibble tibble as_tibble
#' @importFrom purrr imap_dfr
extract_cbps_results <- function(cbps_results,
                                 include_balance_details = TRUE,
                                 smd_threshold = 0.1,
                                 corr_threshold = 0.1,
                                 vr_range = c(0.5, 2),
                                 robust_type = "HC0") {

  required_packages <- c("dplyr", "tibble", "purrr", "sandwich")
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  }

  final_n <- cbps_results$sample_sizes$final

  # Helpers ---------------------------------------------------------------
  safe_min <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
  safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

  # Robust (sandwich) coefficient table for a single weighted GLM
  tidy_robust <- function(model, level = 0.95) {
    est <- stats::coef(model)
    # glm_weightit already returns robust / M-estimation vcov
    V <- if (inherits(model, "glm_weightit")) {
      stats::vcov(model)
    } else {
      sandwich::vcovHC(model, type = robust_type)
    }
    se <- sqrt(diag(V))[names(est)]
    stat <- est / se
    use_t <- identical(model$family$family, "gaussian")
    df_res <- model$df.residual
    alpha <- (1 - level) / 2
    crit <- if (use_t) stats::qt(1 - alpha, df_res) else stats::qnorm(1 - alpha)
    p <- if (use_t) {
      2 * stats::pt(abs(stat), df_res, lower.tail = FALSE)
    } else {
      2 * stats::pnorm(-abs(stat))
    }
    tibble::tibble(
      variable      = names(est),
      glm_estimate  = unname(est),
      glm_se        = unname(se),
      glm_statistic = unname(stat),
      glm_p_value   = unname(p),
      glm_conf_low  = unname(est - crit * se),
      glm_conf_high = unname(est + crit * se),
      se_type       = if (inherits(model, "glm_weightit")) "glm_weightit" else paste0("robust_", robust_type)
    )
  }

  # Pull one imputation's bal.tab Balance table as a tibble
  get_balance_df <- function(bal, imp) {
    if (is.null(bal) || !("Balance" %in% names(bal))) return(NULL)
    df <- as.data.frame(bal$Balance)
    df$variable <- rownames(df)
    df$imputation <- imp
    tibble::as_tibble(df)
  }

  # 1. MODEL ESTIMATES (imputation 1, robust SEs) --------------------------
  model_summary <- tidy_robust(cbps_results$model) %>%
    dplyr::mutate(sample_size = final_n)

  # 2. BOOTSTRAP RESULTS ---------------------------------------------------
  if (!is.null(cbps_results$bootstrap_summary)) {
    bootstrap_results <- cbps_results$bootstrap_summary %>%
      dplyr::rename(dplyr::any_of(c(
        variable            = "term",
        bootstrap_estimate  = "estimate",
        bootstrap_se        = "se",
        bootstrap_conf_low  = "ci_lower",
        bootstrap_conf_high = "ci_upper"
      )))
    combined_estimates <- model_summary %>%
      dplyr::left_join(bootstrap_results, by = "variable")
  } else {
    combined_estimates <- model_summary
  }

  if (!include_balance_details) {
    return(combined_estimates %>%
             dplyr::mutate(
               initial_sample_size = cbps_results$sample_sizes$initial,
               final_sample_size   = final_n,
               n_dropped           = cbps_results$sample_sizes$initial - final_n
             ))
  }

  # 3. BALANCE ACROSS ALL IMPUTATIONS --------------------------------------
  imp_list <- cbps_results$imputation_results
  if (is.null(imp_list) || length(imp_list) == 0) {
    imp_list <- list(list(balance = cbps_results$balance))
  }

  balance_long <- purrr::imap_dfr(imp_list, ~ get_balance_df(.x$balance, .y))

  treatment_type <- "unknown"
  balance_statistic <- NA_character_
  threshold <- NA_real_

  if (nrow(balance_long) > 0) {
    balance_long <- balance_long %>%
      dplyr::select(-dplyr::matches("Threshold")) %>%   # drop cobalt's own flags
      dplyr::rename(dplyr::any_of(c(
        type            = "Type",
        diff_unadj      = "Diff.Un",
        diff_adj        = "Diff.Adj",
        diff_target_adj = "Diff.Target.Adj",  # target-mean difference, NOT a balance SMD
        corr_unadj      = "Corr.Un",
        corr_adj        = "Corr.Adj",
        var_ratio_unadj = "V.Ratio.Un",
        var_ratio_adj   = "V.Ratio.Adj"
      )))

    if ("corr_adj" %in% names(balance_long)) {
      treatment_type <- "continuous"
      balance_statistic <- "|r|"
      threshold <- corr_threshold
      balance_long$balance_stat <- abs(balance_long$corr_adj)
      balance_long$balance_stat_unadj <- if ("corr_unadj" %in% names(balance_long)) abs(balance_long$corr_unadj) else NA_real_
    } else if ("diff_adj" %in% names(balance_long)) {
      treatment_type <- "binary"
      balance_statistic <- "|SMD|"
      threshold <- smd_threshold
      balance_long$balance_stat <- abs(balance_long$diff_adj)
      balance_long$balance_stat_unadj <- if ("diff_unadj" %in% names(balance_long)) abs(balance_long$diff_unadj) else NA_real_
    }

    if (!"var_ratio_adj" %in% names(balance_long)) balance_long$var_ratio_adj <- NA_real_
    if (!"type" %in% names(balance_long)) balance_long$type <- NA_character_
  }

  if (treatment_type != "unknown") {

    balance_df <- balance_long %>%
      dplyr::group_by(variable, type) %>%
      dplyr::summarise(
        n_imputations      = dplyr::n(),
        mean_balance_stat  = safe_mean(balance_stat),
        max_balance_stat   = safe_max(balance_stat),
        max_unadj_stat     = safe_max(balance_stat_unadj),
        min_var_ratio      = safe_min(var_ratio_adj),
        max_var_ratio      = safe_max(var_ratio_adj),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        stat_label = dplyr::case_when(
          treatment_type == "continuous" ~ "|r|",
          type == "Binary"               ~ "|raw diff|",  # cobalt default for binary covariates
          TRUE                           ~ "|SMD|"
        ),
        stat_balanced = max_balance_stat < threshold,  # worst case across imputations
        vr_balanced = if (treatment_type == "binary") {
          is.na(min_var_ratio) | (min_var_ratio >= vr_range[1] & max_var_ratio <= vr_range[2])
        } else {
          NA
        },
        balance_status = dplyr::case_when(
          treatment_type == "continuous" &  stat_balanced ~ paste0("Balanced (|r| < ", threshold, ")"),
          treatment_type == "continuous"                  ~ paste0("Imbalanced (|r| >= ", threshold, ")"),
          stat_balanced & vr_balanced                     ~ "Balanced",
          stat_balanced & !vr_balanced                    ~ "SMD balanced; VR outside range",
          TRUE                                            ~ paste0("Imbalanced (|SMD| >= ", threshold, ")")
        ),
        sample_size = final_n
      )

    summary_stats <- tibble::tibble(
      metric = c(
        "n_imputations",
        "n_covariates",
        "n_covariates_balanced",
        "max_abs_balance_stat_after",
        "mean_abs_balance_stat_after",
        "max_abs_balance_stat_before",
        "min_var_ratio_after",
        "max_var_ratio_after"
      ),
      value = c(
        length(unique(balance_long$imputation)),
        nrow(balance_df),
        sum(balance_df$balance_status %in% c("Balanced", paste0("Balanced (|r| < ", threshold, ")"))),
        safe_max(balance_long$balance_stat),
        safe_mean(balance_long$balance_stat),
        safe_max(balance_long$balance_stat_unadj),
        safe_min(balance_long$var_ratio_adj),
        safe_max(balance_long$var_ratio_adj)
      )
    )

  } else {
    balance_df <- tibble::tibble(variable = "Balance data not available", sample_size = final_n)
    summary_stats <- NULL
  }

  # 4. POOLED RESULTS (PRIMARY) ---------------------------------------------
  pooled_estimates <- if (!is.null(cbps_results$pooled_results)) {
    cbps_results$pooled_results %>%
      dplyr::rename(dplyr::any_of(c(
        variable         = "term",
        pooled_estimate  = "estimate",
        pooled_se        = "se",
        pooled_statistic = "statistic",
        pooled_p_value   = "p.value",
        pooled_conf_low  = "ci_lower",
        pooled_conf_high = "ci_upper",
        pooled_df        = "df",
        pooled_fmi       = "fmi"
      )))
  } else {
    NULL
  }

  # 5. RETURN ----------------------------------------------------------------
  list(
    estimates             = combined_estimates,
    balance               = balance_df,
    balance_by_imputation = balance_long,
    sample_sizes = tibble::tibble(
      metric = c("initial_n", "final_n", "n_dropped"),
      value  = c(
        cbps_results$sample_sizes$initial,
        final_n,
        cbps_results$sample_sizes$initial - final_n
      )
    ),
    summary_stats     = summary_stats,
    pooled_estimates  = pooled_estimates,
    treatment_type    = treatment_type,
    balance_statistic = balance_statistic
  )
}


#' Quick Results Summary for Visualization
#'
#' @param cbps_results List. Output from \code{npcbps_weighted_analysis()}
#' @export
quick_results_summary <- function(cbps_results) {
  result <- extract_cbps_results(cbps_results, include_balance_details = FALSE)

  desired_cols <- c(
    "variable",
    "glm_estimate", "glm_se", "glm_p_value", "glm_conf_low", "glm_conf_high", "se_type",
    "bootstrap_estimate", "bootstrap_se", "bootstrap_conf_low", "bootstrap_conf_high",
    "final_sample_size"
  )

  result %>% dplyr::select(dplyr::any_of(desired_cols))
}
