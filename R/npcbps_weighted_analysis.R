#' Nonparametric CBPS Propensity Score Weighted GLM Analysis with Proper Multiple Imputation
#'
#' Performs comprehensive Nonparametric Covariate Balancing Propensity Score (NPCBPS) weighted analysis
#' with proper multiple imputation using Rubin's rules for pooling.
#'
#' @param data A data frame containing the analysis variables
#' @param outcome_var Character. Name of the outcome variable (y-variable)
#' @param treatment_var Character. Name of the treatment variable (x-variable)
#' @param additional_predictors Character vector. Additional fixed effects/random effects to include in the outcome model
#' @param imputation_vars Character vector. Variables to impute missing values for
#' @param imputation_predictors Character vector. Variables to help with imputation but not used in final model
#' @param propensity_covariates Character vector. Variables to include in propensity weighting (in addition to imputed vars)
#' @param mice_m Integer. Number of multiple imputations (default: 5)
#' @param mice_method Character. MICE imputation method (default: "pmm" for Predictive Mean Matching)
#' @param mice_seed Integer. Random seed for MICE imputation (default: 500)
#' @param cbps_estimand Character. CBPS estimand: "ATE" (Average Treatment Effect) or "ATT" (Average Treatment on Treated) (default: "ATE")
#' @param balance_threshold_m Numeric. Balance threshold for mean differences (default: 0.1)
#' @param balance_threshold_v Numeric. Balance threshold for variance ratios (default: 2)
#' @param family A family object specifying the error distribution and link function for GLM (default: gaussian())
#' @param verbose Logical. Whether to print progress messages (default: TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{pooled_results}: Pooled coefficients using Rubin's rules
#'   \item \code{imputation_results}: List of results from each imputation
#'   \item \code{balance_summary}: Summary of balance across imputations
#'   \item \code{sample_sizes}: List with initial and final sample sizes
#' }
#'
#' @details
#' This function implements the correct multiple imputation workflow:
#' \enumerate{
#'   \item Data cleaning and missing data assessment
#'   \item Multiple imputation using MICE with MCAR testing
#'   \item For EACH imputation:
#'     \itemize{
#'       \item Estimate Nonparametric CBPS propensity scores
#'       \item Assess covariate balance
#'       \item Fit weighted GLM
#'     }
#'   \item Pool results across imputations using Rubin's rules
#' }
#'
#' Rubin's rules combine:
#' - Within-imputation variance (average variance across imputations)
#' - Between-imputation variance (variance of estimates across imputations)
#' - Simulation variance due to finite number of imputations
#'
#' @references
#' Rubin, D. B. (1987). Multiple Imputation for Nonresponse in Surveys. Wiley.
#'
#' van Buuren, S. (2018). Flexible Imputation of Missing Data (2nd ed.). CRC Press.
#'
#' @export
npcbps_weighted_analysis <- function(
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
    balance_threshold_m = 0.1,
    balance_threshold_v = 2,
    family = gaussian(),
    verbose = TRUE
) {

  # Load required libraries
  required_packages <- c("dplyr", "mice", "purrr", "WeightIt", "cobalt", "tibble")
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  }

  # Handle NULL propensity_covariates
  if (is.null(propensity_covariates)) {
    propensity_covariates <- character(0)
  }

  if (verbose) cat("Starting Nonparametric CBPS weighted analysis with proper MI...\n")

  # 1. DATA PREPARATION AND CLEANING
  if (verbose) cat("Step 1: Preparing and cleaning data...\n")

  # Function to extract variable names from formula terms
  extract_vars_from_formula <- function(terms) {
    all_vars <- character(0)
    for (term in terms) {
      vars <- unlist(strsplit(term, "[*:]"))
      vars <- trimws(vars)
      all_vars <- c(all_vars, vars)
    }
    return(unique(all_vars))
  }

  # Extract actual variable names from additional_predictors
  if (!is.null(additional_predictors)) {
    actual_predictor_vars <- extract_vars_from_formula(additional_predictors)
  } else {
    actual_predictor_vars <- character(0)
  }

  # Build variable selection list
  all_vars <- c(outcome_var, treatment_var, actual_predictor_vars,
                imputation_vars, imputation_predictors, propensity_covariates)
  all_vars <- unique(all_vars[all_vars != ""])

  # Select and clean data
  df_clean <- data %>%
    dplyr::select(all_of(all_vars))

  if (verbose) cat(paste("Sample size after variable selection:", nrow(df_clean), "\n"))

  # Remove rows where key variables are missing
  df_clean <- df_clean %>%
    filter(
      !is.na(.data[[outcome_var]]),
      !is.na(.data[[treatment_var]])
    )

  if (verbose) cat(paste("Sample size after removing missing outcome/treatment:", nrow(df_clean), "\n"))

  # Remove additional_predictors that are NA
  if (!is.null(additional_predictors)) {
    for (var in actual_predictor_vars) {
      if (var != treatment_var) {
        n_before <- nrow(df_clean)
        df_clean <- df_clean %>% filter(!is.na(.data[[var]]))
        n_after <- nrow(df_clean)
        if (verbose && n_before != n_after) {
          cat(paste("Removed", n_before - n_after, "rows due to missing", var, "\n"))
        }
      }
    }
  }

  initial_n <- nrow(df_clean)
  if (verbose) cat(paste("Initial sample size after cleaning:", initial_n, "\n"))

  # 2. MISSING DATA ANALYSIS (same as before)
  if (verbose) cat("Step 2: Analyzing missing data patterns and mechanisms...\n")

  missing_data_check <- df_clean %>%
    dplyr::select(dplyr::all_of(c(imputation_vars, imputation_predictors)))

  if (requireNamespace("naniar", quietly = TRUE)) {
    missing_summary <- missing_data_check %>%
      naniar::miss_var_summary()

    if (verbose && nrow(missing_summary) > 0) {
      cat("  Missing data summary:\n")
      missing_vars <- missing_summary %>% dplyr::filter(as.numeric(n_miss) > 0)
      if (nrow(missing_vars) > 0) {
        for (i in 1:nrow(missing_vars)) {
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
          cat("  Missing data mechanism test (MCAR):\n")
          cat(paste("    Little's MCAR test p-value:", round(mcar_result$p.value, 4), "\n"))
          if (mcar_result$p.value > 0.05) {
            cat("     Data appears to be Missing Completely At Random (MCAR)\n")
          } else if (mcar_result$p.value > 0.01) {
            cat("    Marginal evidence against MCAR - assuming MAR for MI\n")
          } else {
            warning("Strong evidence against MCAR (p < 0.01) - results assume MAR", call. = FALSE)
          }
        }
      }, error = function(e) {
        if (verbose) cat("    Could not perform MCAR test\n")
      })
    }
  }

  # 3. MULTIPLE IMPUTATION
  if (verbose) cat("Step 3: Performing multiple imputation (m = ", mice_m, ")...\n")

  # Prepare imputation dataset
  imputation_dataset <- df_clean %>%
    dplyr::select(all_of(c(imputation_vars, imputation_predictors)))

  # Check if imputation is needed
  imputation_needed <- any(is.na(imputation_dataset %>% dplyr::select(all_of(imputation_vars))))

  if (imputation_needed) {
    set.seed(mice_seed)
    imputed_data <- mice(imputation_dataset,
                         m = mice_m,
                         method = mice_method,
                         seed = mice_seed,
                         printFlag = FALSE)

    # Extract each completed dataset (KEEP THEM SEPARATE - don't average!)
    imputed_datasets <- map(1:mice_m, ~complete(imputed_data, .x) %>%
                             dplyr::select(all_of(imputation_vars)))

    if (verbose) cat("  Created", mice_m, "imputed datasets\n")
  } else {
    # No imputation needed - replicate the same dataset m times
    single_dataset <- imputation_dataset %>% dplyr::select(all_of(imputation_vars))
    imputed_datasets <- map(1:mice_m, ~single_dataset)
    if (verbose) cat("  No missing data - using original dataset\n")
  }

  # 4. ANALYZE EACH IMPUTATION SEPARATELY
  if (verbose) cat("Step 4: Analyzing each imputed dataset separately...\n")

  # Get non-imputed variables
  outcome_treatment_data <- df_clean %>%
    dplyr::select(all_of(c(outcome_var, treatment_var, actual_predictor_vars)))

  # Build formulas
  ps_covariate_vars <- c(imputation_vars, propensity_covariates)
  ps_formula <- as.formula(paste(treatment_var, "~", paste(ps_covariate_vars, collapse = " + ")))

  predictor_vars <- c(treatment_var, additional_predictors)
  outcome_formula <- as.formula(paste(outcome_var, "~", paste(predictor_vars, collapse = " + ")))

  # Storage for results from each imputation
  imputation_results <- list()

  for (imp in 1:mice_m) {
    if (verbose) cat(paste("  Processing imputation", imp, "of", mice_m, "...\n"))

    # Combine non-imputed data with this imputation
    df_imp <- bind_cols(outcome_treatment_data, imputed_datasets[[imp]])

    # Check for missing values in propensity covariates
    missing_check <- df_imp %>%
      dplyr::select(all_of(ps_covariate_vars)) %>%
      summarise(across(everything(), ~sum(is.na(.x))))

    if (any(missing_check > 0)) {
      df_imp <- df_imp %>%
        filter(if_all(all_of(ps_covariate_vars), ~!is.na(.x)))
      if (verbose) cat(paste("    Removed", initial_n - nrow(df_imp), "rows with missing covariates\n"))
    }

    # Fit NPCBPS
    tryCatch({
      npcbps_weights <- WeightIt::weightit(
        ps_formula,
        data = df_imp,
        method = "npcbps",
        estimand = cbps_estimand
      )

      # Balance assessment
      balance_table <- bal.tab(npcbps_weights,
                               treat = treatment_var,
                               method = "weighting",
                               m.threshold = balance_threshold_m,
                               v.threshold = balance_threshold_v)

      # Fit weighted GLM
      weighted_model <- glm(outcome_formula,
                           data = df_imp,
                           weights = npcbps_weights$weights,
                           family = family)

      # Store results
      imputation_results[[imp]] <- list(
        data = df_imp,
        weights = npcbps_weights,
        balance = balance_table,
        model = weighted_model,
        coefficients = coef(weighted_model),
        vcov = vcov(weighted_model),
        n = nrow(df_imp)
      )

      if (verbose) cat(paste("    n =", nrow(df_imp), "- Model converged:", weighted_model$converged, "\n"))

    }, error = function(e) {
      warning(paste("Error in imputation", imp, ":", e$message), call. = FALSE)
      imputation_results[[imp]] <- NULL
    })
  }

  # Remove failed imputations
  imputation_results <- imputation_results[!sapply(imputation_results, is.null)]
  n_successful <- length(imputation_results)

  if (n_successful == 0) {
    stop("All imputations failed. Cannot proceed with analysis.")
  }

  if (n_successful < mice_m) {
    warning(paste("Only", n_successful, "of", mice_m, "imputations succeeded"), call. = FALSE)
  }

  # 5. POOL RESULTS USING RUBIN'S RULES
  if (verbose) cat("Step 5: Pooling results using Rubin's rules...\n")

  # Get coefficient names (excluding intercept for reporting)
  all_coef_names <- names(imputation_results[[1]]$coefficients)
  coef_names_no_intercept <- all_coef_names[all_coef_names != "(Intercept)"]

  # Extract coefficients and variances from each imputation
  Q_m <- sapply(imputation_results, function(x) x$coefficients)  # Coefficient estimates
  U_m <- sapply(imputation_results, function(x) diag(x$vcov))    # Variances

  # Rubin's rules formulas
  # Q_bar: pooled coefficient (average across imputations)
  Q_bar <- rowMeans(Q_m)

  # U_bar: within-imputation variance (average of variances)
  U_bar <- rowMeans(U_m)

  # B: between-imputation variance
  B <- apply(Q_m, 1, var)

  # T: total variance
  # T = U_bar + B + B/m (where the B/m term accounts for finite m)
  m <- n_successful
  T <- U_bar + B + B/m

  # Standard errors
  SE <- sqrt(T)

  # Degrees of freedom (Barnard & Rubin 1999 small-sample adjustment)
  lambda <- (B + B/m) / T  # Fraction of missing information
  df_old <- (m - 1) / lambda^2

  # Observed data degrees of freedom (use first imputation's residual df as approximation)
  df_obs <- imputation_results[[1]]$model$df.residual
  df_adj <- (df_obs + 1) / (df_obs + 3) * df_obs * (1 - lambda)
  df <- (df_old * df_adj) / (df_old + df_adj)

  # Confidence intervals (using t-distribution with adjusted df)
  ci_lower <- Q_bar - qt(0.975, df) * SE
  ci_upper <- Q_bar + qt(0.975, df) * SE

  # P-values (two-tailed)
  t_stat <- Q_bar / SE
  p_value <- 2 * pt(abs(t_stat), df, lower.tail = FALSE)

  # Fraction of Missing Information
  r <- (B + B/m) / U_bar  # Relative increase in variance
  FMI <- (r + 2/(df + 3)) / (r + 1)

  # Create pooled results table
  pooled_results <- tibble(
    term = names(Q_bar),
    estimate = Q_bar,
    se = SE,
    statistic = t_stat,
    p.value = p_value,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    df = df,
    fmi = FMI,
    within_var = U_bar,
    between_var = B,
    total_var = T
  )

  if (verbose) {
    cat("\nPooled Results (Rubin's Rules):\n")
    cat("Number of imputations:", m, "\n")
    cat("\nCoefficients:\n")
    print(pooled_results %>%
            dplyr::select(term, estimate, se, ci_lower, ci_upper, p.value, fmi) %>%
            dplyr::filter(term != "(Intercept)"),
          digits = 4)
  }

  # 6. BALANCE SUMMARY ACROSS IMPUTATIONS
  if (verbose) cat("\nStep 6: Summarizing balance across imputations...\n")

  balance_summary <- map_dfr(1:length(imputation_results), function(i) {
    bal <- imputation_results[[i]]$balance
    if ("Balance" %in% names(bal)) {
      balance_df <- as.data.frame(bal$Balance)

      # Find standardized difference column
      diff_col <- NULL
      if ("Diff.Adj" %in% names(balance_df)) {
        diff_col <- "Diff.Adj"
      } else if ("Diff.Target.Adj" %in% names(balance_df)) {
        diff_col <- "Diff.Target.Adj"
      }

      if (!is.null(diff_col)) {
        tibble(
          imputation = i,
          max_abs_std_diff = max(abs(balance_df[[diff_col]]), na.rm = TRUE),
          mean_abs_std_diff = mean(abs(balance_df[[diff_col]]), na.rm = TRUE),
          n_above_threshold = sum(abs(balance_df[[diff_col]]) > balance_threshold_m, na.rm = TRUE)
        )
      } else {
        NULL
      }
    } else {
      NULL
    }
  })

  if (!is.null(balance_summary) && nrow(balance_summary) > 0) {
    if (verbose) {
      cat("Balance summary across imputations:\n")
      cat(paste("  Average max abs std diff:", round(mean(balance_summary$max_abs_std_diff), 3), "\n"))
      cat(paste("  Average mean abs std diff:", round(mean(balance_summary$mean_abs_std_diff), 3), "\n"))
    }

    # Check if any imputation has poor balance
    poor_balance <- any(balance_summary$max_abs_std_diff > balance_threshold_m)
    if (poor_balance) {
      warning("Poor balance detected in one or more imputations", call. = FALSE)
    }
  }

  # 7. RETURN RESULTS
  final_n <- median(sapply(imputation_results, function(x) x$n))

  if (verbose) cat(paste("\nAnalysis completed. Median sample size across imputations:", final_n, "\n"))

  return(list(
    pooled_results = pooled_results,
    imputation_results = imputation_results,
    balance_summary = balance_summary,
    sample_sizes = list(
      initial = initial_n,
      final_median = final_n,
      final_range = range(sapply(imputation_results, function(x) x$n))
    ),
    n_imputations = m,
    mice_seed = mice_seed
  ))
}
