#' CBPS Propensity Score Weighted GLM Analysis
#'
#' Performs comprehensive Covariate Balancing Propensity Score (CBPS) weighted analysis
#' with multiple imputation, balance assessment, and bootstrap confidence intervals.
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
#' @param cbps_stop_method Character. Balance criteria for tuning (default: "es.mean")
#' @param bootstrap_n Integer. Number of bootstrap samples (default: 1000)
#' @param bootstrap_seed Integer. Random seed for bootstrap (default: 20250417)
#' @param balance_threshold_m Numeric. Balance threshold for mean differences (default: 0.1)
#' @param balance_threshold_v Numeric. Balance threshold for variance ratios (default: 2)
#' @param verbose Logical. Whether to print progress messages (default: TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{data}: Final analysis dataset with weights
#'   \item \code{model}: Fitted weighted GLM model
#'   \item \code{weights}: CBPS weight object from WeightIt
#'   \item \code{balance}: Balance assessment from cobalt
#'   \item \code{bootstrap_summary}: Bootstrap confidence intervals
#'   \item \code{bootstrap_results}: Full bootstrap results matrix
#'   \item \code{sample_sizes}: List with initial and final sample sizes
#' }
#'
#' @details
#' This function implements a comprehensive workflow for CBPS analysis including:
#' \enumerate{
#'   \item Data cleaning and missing data assessment
#'   \item Multiple imputation using MICE with MCAR testing
#'   \item CBPS propensity score estimation
#'   \item Covariate balance assessment
#'   \item Weighted GLM fitting
#'   \item Bootstrap confidence intervals
#' }
#'
#' The function tests for Missing Completely At Random (MCAR) assumptions and provides
#' warnings about potential Missing At Random (MAR) or Missing Not At Random (MNAR) scenarios.
#'
#' @references
#' Allison, P. (2015). Imputation by predictive mean matching: Promise & peril. Statistical Horizons.
#'
#' Austin, P. C. (2009). Balance diagnostics for comparing the distribution of baseline
#' covariates between treatment groups in propensity-score matched samples. Statistics in Medicine, 28(25), 3083-107.
#'
#' Li, Y., & Li, L. (2021). Propensity score analysis methods with balancing constraints:
#' A Monte Carlo study. Statistical Methods in Medical Research, 30(4), 1119-1142.
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' results <- cbps_weighted_analysis(
#'   data = my_data,
#'   outcome_var = "weight_percentile",
#'   treatment_var = "treatment_group",
#'   additional_predictors = c("age", "gender"),
#'   imputation_vars = c("mother_age", "gestational_age"),
#'   verbose = TRUE
#' )
#'
#' # Access results
#' summary(results$model)
#' results$bootstrap_summary
#' results$balance
#' }
#'
#' @export
#' @importFrom dplyr select all_of filter if_all mutate bind_cols group_by summarise across rowwise ungroup row_number left_join
#' @importFrom mice mice complete
#' @importFrom purrr map
#' @importFrom WeightIt weightit
#' @importFrom cobalt bal.tab
#' @importFrom tibble tibble
#' @importFrom stats glm as.formula coef quantile sd
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
    cbps_stop_method = "es.mean",
    bootstrap_n = 1000,
    bootstrap_seed = 20250417,
    balance_threshold_m = 0.1,
    balance_threshold_v = 2,
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

  if (verbose) cat("Starting CBPS weighted analysis...\n")

  # 1. DATA PREPARATION AND CLEANING
  if (verbose) cat("Step 1: Preparing and cleaning data...\n")

  # Build variable selection list
  all_vars <- c(outcome_var, treatment_var, additional_predictors,
                imputation_vars, imputation_predictors, propensity_covariates)

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

  # Remove additional_predictors that are NA if specified
  if (!is.null(additional_predictors)) {
    for (var in additional_predictors) {
      n_before <- nrow(df_clean)
      df_clean <- df_clean %>% filter(!is.na(.data[[var]]))
      n_after <- nrow(df_clean)
      if (verbose && n_before != n_after) {
        cat(paste("Removed", n_before - n_after, "rows due to missing", var, "\n"))
      }
    }
  }

  initial_n <- nrow(df_clean)
  if (verbose) cat(paste("Initial sample size after cleaning:", initial_n, "\n"))

  # 2. MISSING DATA ANALYSIS
  if (verbose) cat("Step 2: Analyzing missing data patterns and mechanisms...\n")

  # Prepare dataset for missing data analysis
  missing_data_check <- df_clean %>%
    dplyr::select(dplyr::all_of(c(imputation_vars, imputation_predictors)))

  # Load naniar if available for missing data analysis
  if (requireNamespace("naniar", quietly = TRUE)) {

    # Check overall missingness
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

    # Test for MCAR (Missing Completely At Random)
    if (any(is.na(missing_data_check))) {
      tryCatch({
        mcar_result <- naniar::mcar_test(missing_data_check)

        if (verbose) {
          cat("  Missing data mechanism test (MCAR):\n")
          cat(paste("    Little's MCAR test p-value:", round(mcar_result$p.value, 4), "\n"))

          if (mcar_result$p.value > 0.05) {
            cat("     Data appears to be Missing Completely At Random (MCAR)\n")
            cat("     Multiple imputation assumptions are well supported\n")
          } else if (mcar_result$p.value > 0.01) {
            cat("    Marginal evidence against MCAR (p < 0.05 but > 0.01)\n")
            cat("    Multiple imputation may still be appropriate, but consider MAR assumption\n")
          } else {
            warning("MISSING DATA WARNING: Strong evidence against MCAR assumption (p < 0.01)\n",
                    "  Data may be Missing At Random (MAR) or Missing Not At Random (MNAR)\n",
                    "   Multiple imputation assumes MAR - results may be biased if MNAR\n",
                    "  Consider sensitivity analyses or alternative approaches",
                    call. = FALSE)
          }
        }

        # Additional MAR vs MNAR assessment
        if (verbose) {
          cat("  Assessing MAR vs MNAR likelihood:\n")

          # Check for extreme missingness percentages (potential MNAR indicator)
          high_missing_vars <- missing_summary %>%
            dplyr::filter(as.numeric(n_miss) > 0, as.numeric(pct_miss) > 50)

          if (nrow(high_missing_vars) > 0) {
            warning("POTENTIAL MNAR WARNING: Variables with >50% missing data detected:\n",
                    paste("  ", as.character(high_missing_vars$variable), " (",
                          round(as.numeric(high_missing_vars$pct_miss), 1), "% missing)", collapse = "\n"),
                    "\n   High missingness may indicate Missing Not At Random (MNAR)",
                    "\n   Consider whether missingness is related to unobserved values",
                    "\n   Examples: income (high earners don't report), sensitive topics, etc.",
                    call. = FALSE)
          }

          # Domain-specific MNAR warnings based on variable types
          potentially_mnar_vars <- missing_summary %>%
            dplyr::filter(as.numeric(n_miss) > 0) %>%
            dplyr::filter(
              grepl("income|salary|wage|earn", tolower(variable)) |
                grepl("weight|bmi|height", tolower(variable)) |
                grepl("age|birth", tolower(variable)) |
                grepl("alcohol|smoke|drug", tolower(variable)) |
                grepl("mental|depression|anxiety", tolower(variable))
            )

          if (nrow(potentially_mnar_vars) > 0) {
            cat("     Variables with potential MNAR risk detected:\n")
            for (i in 1:nrow(potentially_mnar_vars)) {
              var_name <- as.character(potentially_mnar_vars$variable[i])
              var_pct <- as.numeric(potentially_mnar_vars$pct_miss[i])

              mnar_reason <- dplyr::case_when(
                grepl("income|salary|wage|earn", tolower(var_name)) ~ "(high earners may not report)",
                grepl("weight|bmi", tolower(var_name)) ~ "(individuals may not report high weights)",
                grepl("age", tolower(var_name)) ~ "(older individuals may not report age)",
                grepl("alcohol|smoke|drug", tolower(var_name)) ~ "(social desirability bias)",
                grepl("mental|depression|anxiety", tolower(var_name)) ~ "(stigma-related non-response)",
                TRUE ~ "(domain-specific non-response patterns)"
              )

              cat(paste("      ", var_name, " (", round(var_pct, 1), "% missing)", mnar_reason, "\n"))
            }
            cat("     Consider sensitivity analyses or domain expert consultation\n")
            cat("     Alternative approaches: selection models, pattern-mixture models\n")
          }

          # Overall MAR vs MNAR assessment
          total_missing_pct <- mean(as.numeric(missing_summary$pct_miss))
          if (total_missing_pct > 20) {
            cat("     Overall high missingness detected (", round(total_missing_pct, 1), "% average)\n")
            cat("     Increased risk of MNAR mechanisms\n")
            cat("     Strongly recommend sensitivity analyses\n")
          } else if (mcar_result$p.value <= 0.05) {
            cat("     Data likely Missing At Random (MAR) given MCAR test results\n")
            cat("     Multiple imputation assumptions are reasonable\n")
          } else {
            cat("     Data appears consistent with MAR assumptions\n")
            cat("    Multiple imputation is well-justified\n")
          }
        }

      }, error = function(e) {
        if (verbose) cat("    Could not perform MCAR test:", e$message, "\n")
      })

    } else {
      if (verbose) cat("    No missing data detected in analysis variables\n")
    }

  } else {
    if (verbose) cat("    naniar package not available - skipping detailed missing data analysis\n")
  }

  # 3. IMPUTATION
  if (verbose) cat("Step 3: Performing multiple imputation...\n")

  # Prepare imputation dataset - include imputation vars AND predictors for better imputation
  imputation_dataset <- df_clean %>%
    dplyr::select(all_of(c(imputation_vars, imputation_predictors))) %>%
    dplyr::mutate(row_id = row_number())

  # Check if imputation is needed
  if (any(is.na(imputation_dataset %>% dplyr::select(all_of(imputation_vars))))) {
    # Perform MICE imputation using both target vars and predictors
    set.seed(mice_seed)
    imputed_data <- mice(imputation_dataset %>% dplyr::select(-row_id),
                         m = mice_m, method = mice_method,
                         seed = mice_seed, printFlag = FALSE)

    # Process imputed data - only extract the imputation target variables
    completed_list <- map(1:mice_m, ~complete(imputed_data, .x) %>%
                            mutate(.imp = .x))

    # Stack and average imputations - only for the target imputation variables
    long_imputed <- bind_rows(completed_list) %>%
      mutate(row_id = rep(1:nrow(imputation_dataset), mice_m)) %>%
      dplyr::select(.imp, row_id, all_of(imputation_vars))

    # Average across imputations for target variables only
    avg_imputations <- long_imputed %>%
      group_by(row_id) %>%
      summarise(
        across(all_of(imputation_vars), ~ mean(.x, na.rm = TRUE),
               .names = "{.col}_imp"),
        .groups = "drop"
      )

    # Join back and replace missing values in target variables only
    df_imp_final <- imputation_dataset %>%
      left_join(avg_imputations, by = "row_id") %>%
      rowwise() %>%
      mutate(
        across(all_of(imputation_vars), ~ ifelse(is.na(.x),
                                                 get(paste0(cur_column(), "_imp")), .x))
      ) %>%
      ungroup() %>%
      dplyr::select(all_of(imputation_vars))  # Only keep the imputed target variables

    if (verbose) cat("Imputation completed.\n")
  } else {
    df_imp_final <- imputation_dataset %>% dplyr::select(all_of(imputation_vars))
    if (verbose) cat("No missing data found in target variables - skipping imputation.\n")
  }

  # 3. COMBINE DATA
  if (verbose) cat("Step 3: Combining imputed and original data...\n")

  outcome_treatment_data <- df_clean %>%
    dplyr::select(all_of(c(outcome_var, treatment_var, additional_predictors)))

  df_final <- bind_cols(outcome_treatment_data, df_imp_final)

  # Build formula for propensity score model - only imputed vars + any additional propensity covariates
  ps_covariate_vars <- c(imputation_vars, propensity_covariates)
  ps_formula <- as.formula(paste(treatment_var, "~", paste(ps_covariate_vars, collapse = " + ")))

  # Check for any remaining missing values in final covariates for propensity model
  missing_check <- df_final %>%
    dplyr::select(all_of(ps_covariate_vars)) %>%
    summarise(across(everything(), ~sum(is.na(.x))))

  if (any(missing_check > 0)) {
    if (verbose) {
      cat("Warning: Missing values detected in covariates after imputation:\n")
      missing_vars <- names(missing_check)[missing_check > 0]
      for (var in missing_vars) {
        cat(paste("  ", var, ":", missing_check[[var]], "missing values\n"))
      }
      cat("Removing rows with missing covariate values...\n")
    }

    # Remove rows with any missing covariates
    df_final <- df_final %>%
      filter(if_all(all_of(ps_covariate_vars), ~!is.na(.x)))

    if (verbose) cat(paste("Sample size after removing missing covariates:", nrow(df_final), "\n"))
  }

  # 5. PROPENSITY SCORE WEIGHTING
  if (verbose) cat("Step 5: Computing CBPS propensity scores...\n")

  cbps_weights <- WeightIt::weightit(
    ps_formula,
    data = df_final,
    method = "cbps",
    estimand = cbps_estimand,
    stop.method = cbps_stop_method
  )

  # 6. BALANCE ASSESSMENT
  if (verbose) cat("Step 6: Assessing covariate balance...\n")

  balance_table <- bal.tab(cbps_weights, treat = treatment_var, method = "weighting",
                           m.threshold = balance_threshold_m, v.threshold = balance_threshold_v)

  # Check balance and provide warnings
  if ("Balance" %in% names(balance_table)) {
    balance_df <- as.data.frame(balance_table$Balance)

    # Check for standardized mean differences - try different column names
    diff_col <- NULL
    if ("Diff.Adj" %in% names(balance_df)) {
      diff_col <- "Diff.Adj"
    } else if ("Diff.Target.Adj" %in% names(balance_df)) {
      diff_col <- "Diff.Target.Adj"
    }

    if (!is.null(diff_col)) {
      max_std_diff <- max(abs(balance_df[[diff_col]]), na.rm = TRUE)
      mean_std_diff <- mean(abs(balance_df[[diff_col]]), na.rm = TRUE)
      unbalanced_vars <- rownames(balance_df)[abs(balance_df[[diff_col]]) > balance_threshold_m]

      if (verbose) {
        cat("Balance Assessment:\n")
        cat(paste("  Maximum absolute standardized difference:", round(max_std_diff, 3), "\n"))
        cat(paste("  Mean absolute standardized difference:", round(mean_std_diff, 3), "\n"))
        cat(paste("  Variables above threshold:", length(unbalanced_vars), "out of", nrow(balance_df), "\n"))
      }

      # Warning for poor balance
      if (max_std_diff > balance_threshold_m) {
        warning(paste("WARNING: Poor covariate balance detected!",
                      "\n  Maximum standardized difference:", round(max_std_diff, 3),
                      "(threshold:", balance_threshold_m, ")",
                      "\n  Variables with poor balance:", paste(unbalanced_vars, collapse = ", "),
                      "\n  Consider different propensity score method or additional covariates."),
                call. = FALSE)
      } else {
        if (verbose) cat("Good covariate balance achieved (all variables < threshold)\n")
      }

      # Additional warning for very poor balance
      if (max_std_diff > 0.25) {
        warning(paste("SEVERE WARNING: Very poor covariate balance detected!",
                      "\n  Maximum standardized difference:", round(max_std_diff, 3),
                      "\n  Results may be unreliable. Consider:",
                      "\n  - Adding more covariates to propensity model",
                      "\n  - Using different weighting method",
                      "\n  - Checking for overlap in propensity scores"),
                call. = FALSE)
      }
    } else {
      if (verbose) cat("Could not find standardized difference column in balance table\n")
    }

    # Check variance ratios if available
    if ("V.Ratio.Adj" %in% names(balance_df)) {
      extreme_var_ratios <- rownames(balance_df)[balance_df$V.Ratio.Adj > balance_threshold_v |
                                                   balance_df$V.Ratio.Adj < (1/balance_threshold_v)]

      if (length(extreme_var_ratios) > 0) {
        warning(paste("WARNING: Extreme variance ratios detected for variables:",
                      paste(extreme_var_ratios, collapse = ", "),
                      "\n  This suggests poor balance in variable distributions."),
                call. = FALSE)
      }
    }
  } else {
    if (verbose) cat("Balance table structure not recognized - skipping balance warnings\n")
  }

  # 7. WEIGHTED GLM
  if (verbose) cat("Step 7: Fitting weighted GLM...\n")

  # Build outcome model formula
  predictor_vars <- c(treatment_var, additional_predictors)
  outcome_formula <- as.formula(paste(outcome_var, "~", paste(predictor_vars, collapse = " + ")))

  weighted_model <- glm(outcome_formula,
                        data = df_final,
                        weights = cbps_weights$weights)

  # 8. BOOTSTRAPPING (Robust version)
  if (verbose) cat("Step 8: Performing bootstrap analysis...\n")

  boot_data <- df_final %>%
    mutate(weight = cbps_weights$weights)

  set.seed(bootstrap_seed)

  # Robust bootstrap function with error handling
  boot_fun <- function(data, indices) {
    tryCatch({
      d <- data[indices, ]

      # Check if bootstrap sample has sufficient variation
      if (length(unique(d[[treatment_var]])) < 2) {
        return(rep(NA, length(predictor_vars)))
      }

      model <- glm(outcome_formula, data = d, weights = weight)

      # Check if model converged
      if (!model$converged) {
        return(rep(NA, length(predictor_vars)))
      }

      coeffs <- coef(model)
      result <- coeffs[predictor_vars]

      # Ensure we return the right length vector
      if (length(result) != length(predictor_vars)) {
        return(rep(NA, length(predictor_vars)))
      }

      return(as.numeric(result))

    }, error = function(e) {
      return(rep(NA, length(predictor_vars)))
    })
  }

  # Run bootstrap with better error handling
  boot_results <- replicate(bootstrap_n, {
    sample_idx <- sample(1:nrow(boot_data), replace = TRUE)
    boot_fun(boot_data, sample_idx)
  }, simplify = FALSE)  # Don't simplify initially

  # Convert to matrix and handle failures
  boot_matrix <- do.call(cbind, boot_results)
  rownames(boot_matrix) <- predictor_vars

  # Check if we have any successful bootstrap samples
  successful_boots <- apply(boot_matrix, 2, function(x) !all(is.na(x)))
  n_successful <- sum(successful_boots)

  if (verbose) cat(paste("Successful bootstrap samples:", n_successful, "out of", bootstrap_n, "\n"))

  if (n_successful < 10) {
    warning("Very few successful bootstrap samples (", n_successful, "). Results may be unreliable.")
  }

  # Calculate bootstrap summary only from successful samples
  if (n_successful > 0) {
    boot_summary <- tibble(
      term = predictor_vars,
      estimate = rowMeans(boot_matrix[, successful_boots, drop = FALSE], na.rm = TRUE),
      se = apply(boot_matrix[, successful_boots, drop = FALSE], 1, sd, na.rm = TRUE),
      ci_lower = apply(boot_matrix[, successful_boots, drop = FALSE], 1, quantile, probs = 0.025, na.rm = TRUE),
      ci_upper = apply(boot_matrix[, successful_boots, drop = FALSE], 1, quantile, probs = 0.975, na.rm = TRUE)
    )
  } else {
    # Fallback if no bootstrap samples succeeded
    warning("No successful bootstrap samples. Using GLM standard errors.")
    boot_summary <- tibble(
      term = predictor_vars,
      estimate = coef(weighted_model)[predictor_vars],
      se = summary(weighted_model)$coefficients[predictor_vars, "Std. Error"],
      ci_lower = NA_real_,
      ci_upper = NA_real_
    )
  }

  final_n <- nrow(df_final)
  if (verbose) cat(paste("Analysis completed. Final sample size:", final_n, "\n"))

  # 9. RETURN RESULTS
  return(list(
    data = df_final,
    model = weighted_model,
    weights = cbps_weights,
    balance = balance_table,
    bootstrap_summary = boot_summary,
    bootstrap_results = boot_matrix,
    sample_sizes = list(initial = initial_n, final = final_n)
  ))
}  # This is the main function closing brace
