#' Extract NPCBPS Results into Structured Format
#'
#' Extracts key results from \code{npcbps_weighted_analysis()} into a structured
#' format suitable for tables, visualization, and further analysis.
#'
#' REVISION NOTES (changes for the updated npcbps_weighted_analysis):
#' \itemize{
#'   \item NEW: returns the pooled average marginal effect (\code{pooled_ame}) and the
#'     per-imputation AMEs (\code{ame_by_imputation}). AME columns are also joined onto
#'     the treatment row of \code{estimates}, so \code{quick_results_summary()} shows them.
#'   \item CHANGED: single-imputation SEs now reuse the vcov the analysis function already
#'     computed (\code{imputation_results[[1]]$vcov}). This works for glm families and for
#'     ordinal (polr) models, where \code{sandwich::vcovHC()} is not available. Falls back
#'     to recomputing for older results objects that lack a stored vcov.
#'   \item NEW: returns \code{family_label} and \code{coef_scale} ("response" for gaussian,
#'     "link" otherwise) so tables can say whether coefficients are log-odds.
#'   \item Works with results from the previous version of npcbps_weighted_analysis
#'     (no AME): AME outputs are then NULL / NA.
#' }
#'
#' Balance is summarized across ALL imputations using the statistic appropriate
#' to the treatment type:
#' \itemize{
#'   \item Binary treatment: |SMD| (primary; raw difference for binary covariates)
#'     plus variance ratio (secondary)
#'   \item Continuous treatment: |r|, the weighted treatment-covariate correlation
#' }
#' Pooled (Rubin's rules) estimates remain the primary results; for non-gaussian
#' outcomes, report \code{pooled_ame} as the effect on the outcome scale.
#'
#' @param cbps_results List. Output from \code{npcbps_weighted_analysis()}
#' @param include_balance_details Logical. Return full list (TRUE) or a simple
#'   estimates data frame (FALSE). Default TRUE.
#' @param smd_threshold Numeric. |SMD| threshold for binary treatments (default 0.1)
#' @param corr_threshold Numeric. |r| threshold for continuous treatments (default 0.1)
#' @param vr_range Numeric length 2. Acceptable variance ratio range (default c(0.5, 2))
#' @param robust_type Character. Sandwich estimator type, used only when the results
#'   object has no stored vcov (default "HC0")
#'
#' @return If \code{include_balance_details = TRUE}, a list containing:
#' \itemize{
#'   \item \code{estimates}: Imputation-1 model estimates (robust SEs) + bootstrap,
#'     with pooled AME columns on the treatment row
#'   \item \code{balance}: Per-covariate balance summarized across imputations
#'   \item \code{balance_by_imputation}: One row per covariate per imputation
#'   \item \code{sample_sizes}: Sample size information
#'   \item \code{summary_stats}: Overall balance metrics across imputations
#'   \item \code{pooled_estimates}: Rubin's rules pooled coefficients
#'   \item \code{pooled_ame}: Rubin's rules pooled AME (PRIMARY EFFECT, outcome scale)
#'   \item \code{ame_by_imputation}: AME estimate and variance from each imputation
#'   \item \code{treatment_type}: "binary" or "continuous"
#'   \item \code{balance_statistic}: "|SMD|" or "|r|"
#'   \item \code{family_label}: description of the outcome model
#'   \item \code{coef_scale}: "response" or "link"
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

  # ---- NEW: model description ----
  family_label <- if (!is.null(cbps_results$family_label)) {
    cbps_results$family_label
  } else if (!is.null(cbps_results$model$family)) {
    paste(cbps_results$model$family$family, "with", cbps_results$model$family$link, "link")
  } else {
    NA_character_
  }
  is_gaussian_identity <- !is.na(family_label) && grepl("^gaussian with identity", family_label)
  coef_scale <- if (is_gaussian_identity) "response" else "link"

  # Helpers ---------------------------------------------------------------
  safe_min <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
  safe_max <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

  # ---- CHANGED: coefficient table using the stored vcov when available ----
  tidy_robust <- function(model, V = NULL, se_label = NULL, level = 0.95) {
    est <- stats::coef(model)
    if (is.null(V)) {
      V <- if (inherits(model, "glm_weightit")) {
        stats::vcov(model)
      } else if (inherits(model, "polr")) {
        tryCatch(sandwich::sandwich(model), error = function(e) stats::vcov(model))
      } else {
        sandwich::vcovHC(model, type = robust_type)
      }
      if (is.null(se_label)) {
        se_label <- if (inherits(model, "glm_weightit")) "glm_weightit" else paste0("robust_", robust_type)
      }
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
      se_type       = if (is.null(se_label)) NA_character_ else se_label
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
  stored_vcov <- if (length(cbps_results$imputation_results) > 0) {
    cbps_results$imputation_results[[1]]$vcov
  } else {
    NULL
  }
  stored_se_label <- if (!is.null(cbps_results$pooled_results$se_method)) {
    cbps_results$pooled_results$se_method[1]
  } else {
    NULL
  }

  model_summary <- tidy_robust(cbps_results$model, V = stored_vcov, se_label = stored_se_label) %>%
    dplyr::mutate(sample_size = final_n, coef_scale = coef_scale)

  # 2. BOOTSTRAP RESULTS ---------------------------------------------------
  if (!is.null(cbps_results$bootstrap_summary)) {
    bootstrap_results <- cbps_results$bootstrap_summary %>%
      dplyr::rename(dplyr::any_of(c(
        variable              = "term",
        bootstrap_estimate    = "estimate",
        bootstrap_se          = "se",
        bootstrap_conf_low    = "ci_lower",
        bootstrap_conf_high   = "ci_upper",
        bootstrap_n_successful = "n_successful"
      )))
    combined_estimates <- model_summary %>%
      dplyr::left_join(bootstrap_results, by = "variable")
  } else {
    combined_estimates <- model_summary
  }

  # ---- NEW: pooled AME, joined onto the treatment row ----
  pooled_ame <- NULL
  if (!is.null(cbps_results$pooled_ame)) {
    pooled_ame <- cbps_results$pooled_ame %>%
      dplyr::rename(dplyr::any_of(c(
        ame_term      = "term",
        ame_estimate  = "estimate",
        ame_se        = "se",
        ame_statistic = "statistic",
        ame_p_value   = "p.value",
        ame_conf_low  = "ci_lower",
        ame_conf_high = "ci_upper",
        ame_df        = "df",
        ame_fmi       = "fmi"
      ))) %>%
      dplyr::mutate(variable = sub("^AME_", "", ame_term))

    combined_estimates <- combined_estimates %>%
      dplyr::left_join(
        pooled_ame %>%
          dplyr::select(dplyr::any_of(c("variable", "ame_estimate", "ame_se", "ame_p_value",
                                        "ame_conf_low", "ame_conf_high", "ame_fmi", "ame_type"))),
        by = "variable"
      )
  }

  if (!include_balance_details) {
    return(combined_estimates %>%
             dplyr::mutate(
               family_label        = family_label,
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
      dplyr::select(-dplyr::matches("Threshold")) %>%
      dplyr::rename(dplyr::any_of(c(
        type            = "Type",
        diff_unadj      = "Diff.Un",
        diff_adj        = "Diff.Adj",
        diff_target_adj = "Diff.Target.Adj",
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
          type == "Binary"               ~ "|raw diff|",
          TRUE                           ~ "|SMD|"
        ),
        stat_balanced = max_balance_stat < threshold,
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

  # 4. POOLED RESULTS ---------------------------------------------------------
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
      ))) %>%
      dplyr::mutate(coef_scale = coef_scale)
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
    pooled_ame        = pooled_ame,
    ame_by_imputation = cbps_results$ame_by_imputation,
    treatment_type    = treatment_type,
    balance_statistic = balance_statistic,
    family_label      = family_label,
    coef_scale        = coef_scale
  )
}


#' Quick Results Summary for Visualization
#'
#' @param cbps_results List. Output from \code{npcbps_weighted_analysis()}
#' @export
quick_results_summary <- function(cbps_results) {
  result <- extract_cbps_results(cbps_results, include_balance_details = FALSE)

  desired_cols <- c(
    "variable", "family_label", "coef_scale",
    "glm_estimate", "glm_se", "glm_p_value", "glm_conf_low", "glm_conf_high", "se_type",
    "bootstrap_estimate", "bootstrap_se", "bootstrap_conf_low", "bootstrap_conf_high",
    "bootstrap_n_successful",
    "ame_estimate", "ame_se", "ame_p_value", "ame_conf_low", "ame_conf_high", "ame_type",
    "final_sample_size"
  )

  result %>% dplyr::select(dplyr::any_of(desired_cols))
}
