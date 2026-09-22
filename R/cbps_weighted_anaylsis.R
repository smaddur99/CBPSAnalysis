#' CBPS Propensity Score Weighted GLM Analysis with Proper Multiple Imputation
#'
#' Performs parametric Covariate Balancing Propensity Score (CBPS) weighted analysis
#' with proper multiple imputation (Rubin's rules), balance assessment appropriate to
#' the treatment type, robust standard errors, and a bootstrap of one imputation's GLM.
#'
#' @param data A data frame containing the analysis variables
#' @param outcome_var Character. Name of the outcome variable (y-variable)
#' @param treatment_var Character. Name of the treatment variable (x-variable)
#' @param additional_predictors Character vector. Additional terms for the outcome model
#' @param imputation_vars Character vector. Variables to impute missing values for
#' @param imputation_predictors Character vector. Variables to help with imputation but not used in final model
#' @param propensity_covariates Character vector. Variables to include in propensity weighting (in addition to imputed vars)
#' @param mice_m Integer. Number of multiple imputations (default: 5)
#' @param mice_method Character. MICE imputation method (default: "pmm")
#' @param mice_seed Integer. Random seed for MICE imputation (default: 500)
#' @param cbps_estimand Character. "ATE" or "ATT" (default: "ATE")
#' @param cbps_stop_method Deprecated. Not used by CBPS; retained so existing calls don't break.
#' @param bootstrap_n Integer. Number of bootstrap samples (default: 1000, 0 to skip)
#' @param bootstrap_seed Integer. Random seed for bootstrap (default: 20250417)
#' @param bootstrap_imputation Integer. Which imputation to bootstrap (default: 1)
#' @param balance_threshold_m Numeric. |SMD| threshold, binary treatments (default: 0.1)
#' @param balance_threshold_v Numeric. Variance ratio threshold; VR must fall in [1/v, v] (default: 2)
#' @param balance_threshold_r Numeric. |r| threshold, continuous treatments (default: 0.1)
#' @param robust_type Character. Sandwich SE type (default: "HC0")
#' @param family GLM family (default: gaussian())
#' @param verbose Logical. Print progress messages (default: TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{data}, \code{model}, \code{weights}, \code{balance}: from imputation 1 (compatibility)
#'   \item \code{bootstrap_summary}, \code{bootstrap_results}: bootstrap of the selected imputation's GLM
#'   \item \code{pooled_results}: Rubin's rules pooled coefficients with robust SEs (PRIMARY RESULTS)
#'   \item \code{imputation_results}: per-imputation data, weights, balance, model, robust vcov
#'   \item \code{balance_summary}: balance summary for each imputation
#'   \item \code{treatment_type}: "binary" or "continuous"
#'   \item \code{sample_sizes}: initial and final sample sizes
#' }
#'
#' @details
#' Workflow:
#' \enumerate{
#'   \item Data cleaning and missing-data assessment (Little's MCAR test)
#'   \item Multiple imputation with MICE; the m completed datasets are kept separate
#'   \item For each imputation: estimate CBPS weights, assess balance, fit weighted GLM,
#'     compute robust (sandwich) variance
#'   \item Pool estimates with Rubin's rules (Barnard-Rubin small-sample df)
#'   \item Bootstrap one imputation's GLM to examine the coefficient distribution
#' }
#'
#' Balance statistics depend on treatment type. Binary treatments: |SMD| (primary) and
#' variance ratio (secondary). Continuous treatments: absolute weighted treatment-covariate
#' correlation |r|, because group-based statistics are undefined.
#'
#' @references
#' Austin, P. C. (2009). Balance diagnostics for comparing the distribution of baseline
#' covariates between treatment groups in propensity-score matched samples. Statistics in Medicine, 28(25), 3083-107.
#'
#' Austin, P. C. (2019). Assessing covariate balance when using the generalized propensity
#' score with quantitative or continuous exposures. Statistical Methods in Medical Research, 28(5), 1365-1377.
#'
#' Barnard, J., & Rubin, D. B. (1999). Small-sample degrees of freedom with multiple imputation.
#' Biometrika, 86(4), 948-955.
#'
#' Fong, C., Hazlett, C., & Imai, K. (2018). Covariate balancing propensity score for a
#' continuous treatment. The Annals of Applied Statistics, 12(1), 156-177.
#'
#' Imai, K., & Ratkovic, M. (2014). Covariate balancing propensity score.
#' Journal of the Royal Statistical Society: Series B, 76(1), 243-263.
#'
#' Rubin, D. B. (1987). Multiple Imputation for Nonresponse in Surveys. Wiley.
#'
#' @export
#' @importFrom dplyr select all_of filter if_all mutate bind_cols summarise across everything
#' @importFrom mice mice complete
#' @importFrom purrr map map_dfr
#' @importFrom WeightIt weightit
#' @importFrom cobalt bal.tab
#' @importFrom tibble tibble
#' @importFrom sandwich vcovHC
#' @importFrom stats glm as.formula coef quantile sd vcov qt pt var na.omit
cbps_weighted_analysis <- function(
    data,
    outcome_var,
    treatment_var,
    additional_predictors = NULL,
    imputation_vars = NULL,
    imputation_predictors = NULL,
    propensity_covariates = NULL,
    mice_m = 5,
    mice_method = "pmm",
    mice_seed = 500,
    cbps_estimand = "ATE",
    cbps_stop_method = NULL,
    bootstrap_n = 1000,
    bootstrap_seed = 20250417,
    bootstrap_imputation = 1,
    balance_threshold_m = 0.1,
    balance_threshold_v = 2,
    balance_threshold_r = 0.1,
    robust_type = "HC0",
    family = gaussian(),
    verbose = TRUE
) {

  required_packages <- c("dplyr", "mice", "purrr", "WeightIt", "cobalt", "tibble", "sandwich")
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  }

  if (!is.null(cbps_stop_method)) {
    message("Note: cbps_stop_method is not used by CBPS and has been ignored.")
  }

  if (is.null(propensity_covariates)) propensity_covariates <- character(0)

  if (verbose) cat("Starting CBPS weighted analysis with proper MI...\n")

  # 1. DATA PREPARATION AND CLEANING -------------------------------------------
  if (verbose) cat("Step 1: Preparing and cleaning data...\n")

  extract_vars_from_formula <- function(terms) {
    all_vars <- character(0)
    for (term in terms) {
      vars <- trimws(unlist(strsplit(term, "[*:]")))
      all_vars <- c(all_vars, vars)
    }
    unique(all_vars)
  }

  actual_predictor_vars <- if (!is.null(additional_predictors)) {
    extract_vars_from_formula(additional_predictors)
  } else {
    character(0)
  }

  all_vars <- c(outcome_var, treatment_var, actual_predictor_vars,
                imputation_vars, imputation_predictors, propensity_covariates)
  all_vars <- unique(all_vars[all_vars != ""])

  df_clean <- data %>%
    dplyr::select(all_of(all_vars)) %>%
    filter(!is.na(.data[[outcome_var]]), !is.na(.data[[treatment_var]]))

  if (verbose) cat(paste("Sample size after removing missing outcome/treatment:", nrow(df_clean), "\n"))

  for (var in setdiff(actual_predictor_vars, treatment_var)) {
    n_before <- nrow(df_clean)
    df_clean <- df_clean %>% filter(!is.na(.data[[var]]))
    if (verbose && n_before != nrow(df_clean)) {
      cat(paste("Removed", n_before - nrow(df_clean), "rows due to missing", var, "\n"))
    }
  }

  initial_n <- nrow(df_clean)
  if (verbose) cat(paste("Initial sample size after cleaning:", initial_n, "\n"))

  # Treatment type determines which balance statistic is valid
  n_treat_levels <- length(unique(stats::na.omit(df_clean[[treatment_var]])))
  if (n_treat_levels == 2) {
    treatment_type <- "binary"
  } else if (is.numeric(df_clean[[treatment_var]])) {
    treatment_type <- "continuous"
  } else {
    stop("Treatment must be binary or numeric continuous; multi-category treatments are not supported.")
  }
  if (verbose) cat(paste("Treatment type detected:", treatment_type, "\n"))

  # 2. MISSING DATA ANALYSIS -----------------------------------------------------
  if (verbose) cat("Step 2: Analyzing missing data patterns...\n")

  missing_data_check <- df_clean %>%
    dplyr::select(dplyr::all_of(c(imputation_vars, imputation_predictors)))

  if (requireNamespace("naniar", quietly = TRUE)) {
    missing_summary <- naniar::miss_var_summary(missing_data_check)

    if (verbose && nrow(missing_summary) > 0) {
      missing_vars <- missing_summary %>% dplyr::filter(as.numeric(n_miss) > 0)
      if (nrow(missing_vars) > 0) {
        cat("  Missing data summary:\n")
        for (i in seq_len(nrow(missing_vars))) {
          cat(paste("    ", as.character(missing_vars$variable[i]), ":",
                    as.numeric(missing_vars$n_miss[i]),
                    "missing (", round(as.numeric(missing_vars$pct_miss[i]), 1), "%)\n"))
        }
      } else {
        cat("    No missing data in imputation variables\n")
      }
    }

    if (any(is.na(missing_data_check))) {
      tryCatch({
        mcar_result <- naniar::mcar_test(missing_data_check)
        if (verbose) {
          cat(paste("  Little's MCAR test p-value:", round(mcar_result$p.value, 4), "\n"))
          if (mcar_result$p.value > 0.05) {
            cat("    No evidence against MCAR (note: this does not prove MCAR)\n")
          } else if (mcar_result$p.value > 0.01) {
            cat("    Marginal evidence against MCAR; MI remains appropriate under MAR\n")
          }
        }
        if (mcar_result$p.value <= 0.01) {
          warning("Strong evidence against MCAR (p < 0.01). MI assumes MAR; ",
                  "results may be biased if data are MNAR. Consider sensitivity analyses.",
                  call. = FALSE)
        }
        high_missing_vars <- missing_summary %>%
          dplyr::filter(as.numeric(n_miss) > 0, as.numeric(pct_miss) > 50)
        if (nrow(high_missing_vars) > 0) {
          warning("Variables with >50% missing data: ",
                  paste(as.character(high_missing_vars$variable), collapse = ", "),
                  ". High missingness may indicate MNAR.", call. = FALSE)
        }
      }, error = function(e) {
        if (verbose) cat("    Could not perform MCAR test:", e$message, "\n")
      })
    }
  } else if (verbose) {
    cat("    naniar not available - skipping missing data analysis\n")
  }

  # 3. MULTIPLE IMPUTATION (datasets kept separate - NOT averaged) ---------------
  if (verbose) cat("Step 3: Multiple imputation...\n")

  imputation_dataset <- df_clean %>%
    dplyr::select(all_of(c(imputation_vars, imputation_predictors)))

  imputation_needed <- any(is.na(imputation_dataset %>% dplyr::select(all_of(imputation_vars))))

  if (imputation_needed) {
    set.seed(mice_seed)
    imputed_data <- mice(imputation_dataset, m = mice_m, method = mice_method,
                         seed = mice_seed, printFlag = FALSE)
    imputed_datasets <- map(seq_len(mice_m), ~ complete(imputed_data, .x) %>%
                              dplyr::select(all_of(imputation_vars)))
    if (verbose) cat("  Created", mice_m, "imputed datasets\n")
  } else {
    # No missing data: a single analysis; Rubin's rules reduce to that analysis
    imputed_datasets <- list(imputation_dataset %>% dplyr::select(all_of(imputation_vars)))
    if (verbose) cat("  No missing data - running a single complete-data analysis\n")
  }

  # 4. ANALYZE EACH IMPUTATION ------------------------------------------------------
  if (verbose) cat("Step 4: Analyzing each imputed dataset...\n")

  outcome_treatment_data <- df_clean %>%
    dplyr::select(all_of(unique(c(outcome_var, treatment_var, actual_predictor_vars))))

  ps_covariate_vars <- unique(c(imputation_vars, propensity_covariates))
  ps_formula <- as.formula(paste(treatment_var, "~", paste(ps_covariate_vars, collapse = " + ")))

  outcome_formula <- as.formula(paste(outcome_var, "~",
                                      paste(c(treatment_var, additional_predictors), collapse = " + ")))

  # Balance thresholds and statistic by treatment type
  if (treatment_type == "binary") {
    bal_thresholds <- c(m = balance_threshold_m, v = balance_threshold_v)
    stat_col <- "Diff.Adj"
    stat_threshold <- balance_threshold_m
    stat_label <- "|SMD|"
  } else {
    bal_thresholds <- c(cor = balance_threshold_r)
    stat_col <- "Corr.Adj"
    stat_threshold <- balance_threshold_r
    stat_label <- "|r|"
  }

  check_balance <- function(balance_table, imp) {
    if (!("Balance" %in% names(balance_table))) return(NULL)
    bdf <- as.data.frame(balance_table$Balance)
    if (!(stat_col %in% names(bdf))) {
      warning(paste("Imputation", imp, ": balance column", stat_col, "not found"), call. = FALSE)
      return(NULL)
    }
    stat <- abs(bdf[[stat_col]])
    unbalanced <- rownames(bdf)[!is.na(stat) & stat >= stat_threshold]

    vr_outside <- character(0)
    if (treatment_type == "binary" && "V.Ratio.Adj" %in% names(bdf)) {
      vr <- bdf$V.Ratio.Adj
      vr_outside <- rownames(bdf)[!is.na(vr) & (vr > balance_threshold_v | vr < 1 / balance_threshold_v)]
    }

    if (verbose) {
      cat(paste0("    Balance (", stat_label, "): max = ", round(max(stat, na.rm = TRUE), 4),
                 ", mean = ", round(mean(stat, na.rm = TRUE), 4),
                 ", above threshold = ", length(unbalanced), " of ", nrow(bdf), "\n"))
    }
    if (length(unbalanced) > 0) {
      warning(paste0("Imputation ", imp, ": ", stat_label, " >= ", stat_threshold,
                     " for: ", paste(unbalanced, collapse = ", ")), call. = FALSE)
    }
    if (treatment_type == "binary" && max(stat, na.rm = TRUE) > 0.25) {
      warning(paste("Imputation", imp, ": severe imbalance (|SMD| > 0.25). Results may be unreliable."),
              call. = FALSE)
    }
    if (length(vr_outside) > 0) {
      warning(paste0("Imputation ", imp, ": variance ratio outside [", 1 / balance_threshold_v, ", ",
                     balance_threshold_v, "] for: ", paste(vr_outside, collapse = ", ")), call. = FALSE)
    }

    tibble(
      imputation = imp,
      statistic = stat_label,
      max_abs_stat = max(stat, na.rm = TRUE),
      mean_abs_stat = mean(stat, na.rm = TRUE),
      n_above_threshold = length(unbalanced),
      n_vr_outside = length(vr_outside)
    )
  }

  imputation_results <- vector("list", length(imputed_datasets))

  for (imp in seq_along(imputed_datasets)) {
    if (verbose) cat(paste("  Imputation", imp, "of", length(imputed_datasets), "\n"))

    df_imp <- bind_cols(outcome_treatment_data, imputed_datasets[[imp]])
    df_imp <- df_imp %>% filter(if_all(all_of(ps_covariate_vars), ~ !is.na(.x)))

    imputation_results[[imp]] <- tryCatch({
      cbps_w <- WeightIt::weightit(ps_formula, data = df_imp,
                                   method = "cbps", estimand = cbps_estimand)

      balance_table <- bal.tab(cbps_w, thresholds = bal_thresholds)
      balance_check <- check_balance(balance_table, imp)

      weighted_model <- glm(outcome_formula, data = df_imp,
                            weights = cbps_w$weights, family = family)

      # Robust (sandwich) variance: model-based glm SEs are invalid with PS weights
      robust_vcov <- sandwich::vcovHC(weighted_model, type = robust_type)

      list(
        data = df_imp,
        weights = cbps_w,
        balance = balance_table,
        balance_check = balance_check,
        model = weighted_model,
        coefficients = coef(weighted_model),
        vcov = robust_vcov,
        n = nrow(df_imp)
      )
    }, error = function(e) {
      warning(paste("Error in imputation", imp, ":", e$message), call. = FALSE)
      NULL
    })
  }

  imputation_results <- imputation_results[!vapply(imputation_results, is.null, logical(1))]
  m <- length(imputation_results)
  if (m == 0) stop("All imputations failed. Cannot proceed with analysis.")
  if (imputation_needed && m < mice_m) {
    warning(paste("Only", m, "of", mice_m, "imputations succeeded"), call. = FALSE)
  }

  # 5. POOL WITH RUBIN'S RULES (robust SEs) -------------------------------------------
  if (verbose) cat("Step 5: Pooling results with Rubin's rules...\n")

  coef_names <- names(imputation_results[[1]]$coefficients)
  Q_m <- do.call(cbind, lapply(imputation_results, function(x) x$coefficients[coef_names]))
  U_m <- do.call(cbind, lapply(imputation_results, function(x) diag(x$vcov)[coef_names]))
  rownames(Q_m) <- rownames(U_m) <- coef_names

  Q_bar <- rowMeans(Q_m)
  U_bar <- rowMeans(U_m)
  B <- if (m > 1) apply(Q_m, 1, stats::var) else rep(0, length(Q_bar))
  T_var <- U_bar + (1 + 1 / m) * B
  SE <- sqrt(T_var)

  # Barnard-Rubin degrees of freedom (handles B = 0, e.g., no missing data)
  lambda <- ((1 + 1 / m) * B) / T_var
  df_obs <- imputation_results[[1]]$model$df.residual
  df_adj <- (df_obs + 1) / (df_obs + 3) * df_obs * (1 - lambda)
  df_old <- ifelse(lambda > 0, (m - 1) / lambda^2, Inf)
  df <- ifelse(is.finite(df_old), (df_old * df_adj) / (df_old + df_adj), df_adj)

  t_stat <- Q_bar / SE
  crit <- qt(0.975, df)
  r_ratio <- ((1 + 1 / m) * B) / U_bar
  FMI <- (r_ratio + 2 / (df + 3)) / (r_ratio + 1)

  pooled_results <- tibble(
    term = coef_names,
    estimate = Q_bar,
    se = SE,
    statistic = t_stat,
    p.value = 2 * pt(abs(t_stat), df, lower.tail = FALSE),
    ci_lower = Q_bar - crit * SE,
    ci_upper = Q_bar + crit * SE,
    df = df,
    fmi = FMI,
    within_var = U_bar,
    between_var = B,
    total_var = T_var,
    se_type = paste0("robust_", robust_type)
  )

  if (verbose) {
    cat("\nPooled results (Rubin's rules, robust SEs) - m =", m, "\n")
    print(pooled_results %>%
            dplyr::select(term, estimate, se, ci_lower, ci_upper, p.value, fmi) %>%
            dplyr::filter(term != "(Intercept)"), digits = 4)
  }

  # 6. BOOTSTRAP ONE IMPUTATION'S GLM ----------------------------------------------------
  bootstrap_summary <- NULL
  bootstrap_results <- NULL
  coef_no_int <- coef_names[coef_names != "(Intercept)"]

  if (bootstrap_n > 0) {
    boot_imp <- min(bootstrap_imputation, m)
    if (verbose) cat("\nStep 6: Bootstrapping GLM from imputation", boot_imp, "...\n")

    boot_data <- imputation_results[[boot_imp]]$data %>%
      mutate(weight = imputation_results[[boot_imp]]$weights$weights)

    set.seed(bootstrap_seed)

    boot_fun <- function(d) {
      tryCatch({
        if (length(unique(d[[treatment_var]])) < 2) return(rep(NA_real_, length(coef_no_int)))
        fit <- glm(outcome_formula, data = d, weights = weight, family = family)
        if (!fit$converged) return(rep(NA_real_, length(coef_no_int)))
        as.numeric(coef(fit)[coef_no_int])
      }, error = function(e) rep(NA_real_, length(coef_no_int)))
    }

    boot_list <- replicate(bootstrap_n, {
      boot_fun(boot_data[sample(seq_len(nrow(boot_data)), replace = TRUE), ])
    }, simplify = FALSE)

    boot_matrix <- do.call(cbind, boot_list)
    rownames(boot_matrix) <- coef_no_int
    ok <- apply(boot_matrix, 2, function(x) !all(is.na(x)))

    if (verbose) cat(paste("Successful bootstrap samples:", sum(ok), "out of", bootstrap_n, "\n"))
    if (sum(ok) < 10) warning("Very few successful bootstrap samples (", sum(ok), ").", call. = FALSE)

    if (sum(ok) > 0) {
      bm <- boot_matrix[, ok, drop = FALSE]
      bootstrap_summary <- tibble(
        term = coef_no_int,
        estimate = rowMeans(bm, na.rm = TRUE),
        se = apply(bm, 1, sd, na.rm = TRUE),
        ci_lower = apply(bm, 1, quantile, probs = 0.025, na.rm = TRUE),
        ci_upper = apply(bm, 1, quantile, probs = 0.975, na.rm = TRUE)
      )
    }
    bootstrap_results <- boot_matrix
  }

  # 7. BALANCE SUMMARY ACROSS IMPUTATIONS ----------------------------------------------
  balance_summary <- map_dfr(imputation_results, "balance_check")

  if (verbose && nrow(balance_summary) > 0) {
    cat(paste0("\nBalance across imputations (", stat_label, "): worst max = ",
               round(max(balance_summary$max_abs_stat), 4),
               ", average mean = ", round(mean(balance_summary$mean_abs_stat), 4), "\n"))
  }
  if (nrow(balance_summary) > 0 && any(balance_summary$n_above_threshold > 0)) {
    warning("Imbalance detected in one or more imputations", call. = FALSE)
  }

  final_n <- stats::median(vapply(imputation_results, function(x) x$n, numeric(1)))
  if (verbose) cat(paste("\nAnalysis completed. Median sample size across imputations:", final_n, "\n"))

  # 8. RETURN -------------------------------------------------------------------------------
  list(
    data = imputation_results[[1]]$data,
    model = imputation_results[[1]]$model,
    weights = imputation_results[[1]]$weights,
    balance = imputation_results[[1]]$balance,
    bootstrap_summary = bootstrap_summary,
    bootstrap_results = bootstrap_results,
    sample_sizes = list(initial = initial_n, final = final_n),
    pooled_results = pooled_results,
    imputation_results = imputation_results,
    balance_summary = balance_summary,
    treatment_type = treatment_type,
    n_imputations = m,
    n_bootstraps = bootstrap_n,
    bootstrap_imputation_used = if (bootstrap_n > 0) min(bootstrap_imputation, m) else NA
  )
}
