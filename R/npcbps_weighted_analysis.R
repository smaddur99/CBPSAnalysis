#' Nonparametric CBPS Propensity Score Weighted GLM Analysis with Proper Multiple Imputation
#'
#' Performs comprehensive Nonparametric Covariate Balancing Propensity Score (NPCBPS) weighted analysis
#' with proper multiple imputation using Rubin's rules, plus bootstrap of a single GLM model.
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
#' @param bootstrap_n Integer. Number of bootstrap samples for the GLM (default: 1000, set to 0 to skip)
#' @param bootstrap_seed Integer. Random seed for bootstrap (default: 20250417)
#' @param bootstrap_imputation Integer. Which imputation to use for bootstrap (default: 1 = first imputation)
#' @param balance_threshold_m Numeric. Balance threshold for mean differences (default: 0.1)
#' @param balance_threshold_v Numeric. Balance threshold for variance ratios (default: 2)
#' @param family A family object specifying the error distribution and link function for GLM (default: gaussian())
#' @param verbose Logical. Whether to print progress messages (default: TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{data}: Final analysis dataset from first imputation (for compatibility)
#'   \item \code{model}: Fitted weighted GLM model from first imputation (for compatibility)
#'   \item \code{weights}: NPCBPS weight object from first imputation (for compatibility)
#'   \item \code{balance}: Balance assessment from first imputation (for compatibility)
#'   \item \code{bootstrap_summary}: Bootstrap CIs from the selected imputation's GLM
#'   \item \code{bootstrap_results}: Bootstrap distribution matrix from the selected imputation's GLM
#'   \item \code{pooled_results}: Rubin's rules pooled coefficients (PRIMARY RESULTS)
#'   \item \code{imputation_results}: List of results from each imputation
#'   \item \code{balance_summary}: Summary of balance across imputations
#'   \item \code{sample_sizes}: List with initial and final sample sizes
#' }
#'
#' @details
#' This function implements proper multiple imputation workflow:
#' \enumerate{
#'   \item Data cleaning and missing data assessment with MCAR testing
#'   \item Multiple imputation using MICE (m imputations kept separate)
#'   \item For EACH imputation:
#'     \itemize{
#'       \item Estimate Nonparametric CBPS propensity scores
#'       \item Assess covariate balance
#'       \item Fit weighted GLM
#'     }
#'   \item Pool results using Rubin's rules (PRIMARY inference)
#'   \item Bootstrap ONE of the fitted GLM models to examine distribution
#' }
#'
#' The bootstrap provides the distribution of coefficients for a single imputation's model,
#' while Rubin's rules provides the primary inference accounting for imputation uncertainty.
#'
#' The nonparametric CBPS method uses kernel-based estimation and does not assume
#' a parametric form for the propensity score model. This can be more flexible than
#' parametric CBPS but may require larger sample sizes.
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
#' Barnard, J., & Rubin, D. B. (1999). Small-sample degrees of freedom with multiple imputation.
#' Biometrika, 86(4), 948-955.
#'
#' Fong, C., Hazlett, C., & Imai, K. (2018). Covariate balancing propensity score for a
#' continuous treatment: Application to the efficacy of political advertisements.
#' The Annals of Applied Statistics, 12(1), 156-177.
#'
#' Rubin, D. B. (1987). Multiple Imputation for Nonresponse in Surveys. Wiley.
#'
#' van Buuren, S. (2018). Flexible Imputation of Missing Data (2nd ed.). CRC Press.
#'
#' @examples
#' \dontrun{
#' # Basic usage with continuous outcome (default gaussian family)
#' results <- npcbps_weighted_analysis(
#'   data = my_data,
#'   outcome_var = "weight_percentile",
#'   treatment_var = "treatment_group",
#'   additional_predictors = c("age", "gender"),
#'   imputation_vars = c("mother_age", "gestational_age"),
#'   verbose = TRUE
#' )
#'
#' # Binary outcome (e.g., preterm birth yes/no)
#' results <- npcbps_weighted_analysis(
#'   data = my_data,
#'   outcome_var = "preterm_birth",
#'   treatment_var = "treatment_group",
#'   additional_predictors = c("age", "gender"),
#'   imputation_vars = c("mother_age", "gestational_age"),
#'   family = binomial(),
#'   verbose = TRUE
#' )
#'
#' # With interaction terms
#' results <- npcbps_weighted_analysis(
#'   data = my_data,
#'   outcome_var = "outcome",
#'   treatment_var = "treatment",
#'   additional_predictors = c("var1", "var2", "var1*var2"),
#'   imputation_vars = c("covariate1", "covariate2"),
#'   verbose = TRUE
#' )
#'
#' # Access results
#' results$pooled_results  # PRIMARY - Rubin's rules
#' results$bootstrap_summary  # Bootstrap distribution of one GLM
#' results$balance
#' }
#'
#' @export
#' @importFrom dplyr select all_of filter if_all mutate bind_cols summarise across everything
#' @importFrom mice mice complete
#' @importFrom purrr map map_dfr
#' @importFrom WeightIt weightit
#' @importFrom cobalt bal.tab
#' @importFrom tibble tibble
#' @importFrom stats glm as.formula coef quantile sd vcov qt pt var
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
    bootstrap_n = 1000,
    bootstrap_seed = 20250417,
    bootstrap_imputation = 1,
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
  if (verbose) cat("DEBUG: Using kernel-based nonparametric CBPS estimation\n")

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

  if (verbose) {
    cat("DEBUG: additional_predictors =", paste(additional_predictors, collapse = ", "), "\n")
    cat("DEBUG: actual_predictor_vars =", paste(actual_predictor_vars, collapse = ", "), "\n")
  }

  # Build variable selection list
  all_vars <- c(outcome_var, treatment_var, actual_predictor_vars,
                imputation_vars, imputation_predictors, propensity_covariates)

  if (verbose) cat("DEBUG: all_vars =", paste(all_vars, collapse = ", "), "\n")

  all_vars <- unique(all_vars[all_vars != ""])

  if (verbose) cat("DEBUG: all_vars after unique =", paste(all_vars, collapse = ", "), "\n")

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

  # 2. MISSING DATA ANALYSIS
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

  # 3. MULTIPLE IMPUTATION
  if (verbose) cat("Step 3: Performing multiple imputation (m = ", mice_m, ")...\n")

  imputation_dataset <- df_clean %>%
    dplyr::select(all_of(c(imputation_vars, imputation_predictors)))

  imputation_needed <- any(is.na(imputation_dataset %>% dplyr::select(all_of(imputation_vars))))

  if (imputation_needed) {
    set.seed(mice_seed)
    imputed_data <- mice(imputation_dataset,
                         m = mice_m,
                         method = mice_method,
                         seed = mice_seed,
                         printFlag = FALSE)

    # KEEP IMPUTATIONS SEPARATE - don't average!
    imputed_datasets <- map(1:mice_m, ~complete(imputed_data, .x) %>%
                              dplyr::select(all_of(imputation_vars)))

    if (verbose) cat("  Created", mice_m, "imputed datasets\n")
  } else {
    single_dataset <- imputation_dataset %>% dplyr::select(all_of(imputation_vars))
    imputed_datasets <- map(1:mice_m, ~single_dataset)
    if (verbose) cat("  No missing data - using original dataset\n")
  }

  # 4. ANALYZE EACH IMPUTATION SEPARATELY
  if (verbose) cat("Step 4: Analyzing each imputed dataset separately...\n")

  outcome_treatment_data <- df_clean %>%
    dplyr::select(all_of(c(outcome_var, treatment_var, actual_predictor_vars)))

  ps_covariate_vars <- c(imputation_vars, propensity_covariates)
  ps_formula <- as.formula(paste(treatment_var, "~", paste(ps_covariate_vars, collapse = " + ")))

  predictor_vars <- c(treatment_var, additional_predictors)
  outcome_formula <- as.formula(paste(outcome_var, "~", paste(predictor_vars, collapse = " + ")))

  if (verbose) {
    cat("DEBUG: Outcome formula:", deparse(outcome_formula), "\n")
    cat("DEBUG: GLM family:", family$family, "with", family$link, "link\n")
  }

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
    if (verbose) cat("    Note: NPCBPS uses kernel-based estimation - may take longer than parametric CBPS\n")

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

      # Check balance and provide warnings
      if ("Balance" %in% names(balance_table)) {
        balance_df <- as.data.frame(balance_table$Balance)

        # Check for standardized mean differences
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
            cat("    Balance Assessment:\n")
            cat(paste("      Maximum absolute standardized difference:", round(max_std_diff, 3), "\n"))
            cat(paste("      Mean absolute standardized difference:", round(mean_std_diff, 3), "\n"))
            cat(paste("      Variables above threshold:", length(unbalanced_vars), "out of", nrow(balance_df), "\n"))
          }

          # Warning for poor balance
          if (max_std_diff > balance_threshold_m) {
            warning(paste("WARNING: Poor covariate balance detected in imputation", imp, "!",
                          "\n  Maximum standardized difference:", round(max_std_diff, 3),
                          "(threshold:", balance_threshold_m, ")",
                          "\n  Variables with poor balance:", paste(unbalanced_vars, collapse = ", "),
                          "\n  Consider adding more covariates or checking sample size."),
                    call. = FALSE)
          }

          # Additional warning for very poor balance
          if (max_std_diff > 0.25) {
            warning(paste("SEVERE WARNING: Very poor covariate balance detected in imputation", imp, "!",
                          "\n  Maximum standardized difference:", round(max_std_diff, 3),
                          "\n  Results may be unreliable. Consider:",
                          "\n  - Adding more covariates to propensity model",
                          "\n  - Checking if sample size is sufficient for nonparametric methods",
                          "\n  - Using parametric CBPS instead",
                          "\n  - Checking for overlap in covariate distributions"),
                    call. = FALSE)
          }
        }

        # Check variance ratios if available
        if ("V.Ratio.Adj" %in% names(balance_df)) {
          extreme_var_ratios <- rownames(balance_df)[balance_df$V.Ratio.Adj > balance_threshold_v |
                                                       balance_df$V.Ratio.Adj < (1/balance_threshold_v)]

          if (length(extreme_var_ratios) > 0) {
            warning(paste("WARNING: Extreme variance ratios detected in imputation", imp, "for variables:",
                          paste(extreme_var_ratios, collapse = ", "),
                          "\n  This suggests poor balance in variable distributions."),
                    call. = FALSE)
          }
        }
      }

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

  # Get coefficient names
  all_coef_names <- names(imputation_results[[1]]$coefficients)
  model_coef_names_no_intercept <- all_coef_names[all_coef_names != "(Intercept)"]

  if (verbose) {
    cat("DEBUG: Model coefficient names:", paste(all_coef_names, collapse = ", "), "\n")
    cat("DEBUG: Coefficients for bootstrap:", paste(model_coef_names_no_intercept, collapse = ", "), "\n")
  }

  # 5. POOL RESULTS USING RUBIN'S RULES
  if (verbose) cat("Step 5: Pooling results using Rubin's rules...\n")

  # Extract coefficients and variances from each imputation
  Q_m <- sapply(imputation_results, function(x) x$coefficients)
  U_m <- sapply(imputation_results, function(x) diag(x$vcov))

  # Rubin's rules
  Q_bar <- rowMeans(Q_m)
  U_bar <- rowMeans(U_m)
  B <- apply(Q_m, 1, stats::var)  # FIXED: Explicitly use stats::var

  m <- n_successful
  T <- U_bar + B + B/m
  SE <- sqrt(T)

  # Degrees of freedom
  lambda <- (B + B/m) / T
  df_old <- (m - 1) / lambda^2
  df_obs <- imputation_results[[1]]$model$df.residual
  df_adj <- (df_obs + 1) / (df_obs + 3) * df_obs * (1 - lambda)
  df <- (df_old * df_adj) / (df_old + df_adj)

  # Confidence intervals
  ci_lower <- Q_bar - qt(0.975, df) * SE
  ci_upper <- Q_bar + qt(0.975, df) * SE

  # P-values
  t_stat <- Q_bar / SE
  p_value <- 2 * pt(abs(t_stat), df, lower.tail = FALSE)

  # FMI
  r <- (B + B/m) / U_bar
  FMI <- (r + 2/(df + 3)) / (r + 1)

  # Pooled results table
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
    cat("\nPooled Results (Rubin's Rules - PRIMARY INFERENCE):\n")
    cat("Number of imputations:", m, "\n")
    cat("\nCoefficients:\n")
    print(pooled_results %>%
            dplyr::select(term, estimate, se, ci_lower, ci_upper, p.value, fmi) %>%
            dplyr::filter(term != "(Intercept)"),
          digits = 4)
  }

  # 6. BOOTSTRAP ONE GLM MODEL
  bootstrap_summary <- NULL
  bootstrap_results <- NULL

  if (bootstrap_n > 0) {
    if (verbose) cat("\nStep 6: Bootstrapping GLM from imputation", bootstrap_imputation, "...\n")

    # Use specified imputation for bootstrap
    boot_data <- imputation_results[[bootstrap_imputation]]$data %>%
      mutate(weight = imputation_results[[bootstrap_imputation]]$weights$weights)

    set.seed(bootstrap_seed)

    # Robust bootstrap function with error handling
    boot_fun <- function(data, indices) {
      tryCatch({
        d <- data[indices, ]

        # Check if bootstrap sample has sufficient variation
        if (length(unique(d[[treatment_var]])) < 2) {
          return(rep(NA, length(model_coef_names_no_intercept)))
        }

        model <- glm(outcome_formula, data = d, weights = weight, family = family)

        # Check if model converged
        if (!model$converged) {
          return(rep(NA, length(model_coef_names_no_intercept)))
        }

        coeffs <- coef(model)

        # Extract coefficients using the actual model coefficient names
        result <- coeffs[model_coef_names_no_intercept]

        # Ensure we return the right length vector
        if (length(result) != length(model_coef_names_no_intercept)) {
          return(rep(NA, length(model_coef_names_no_intercept)))
        }

        return(as.numeric(result))

      }, error = function(e) {
        return(rep(NA, length(model_coef_names_no_intercept)))
      })
    }

    # Run bootstrap with better error handling
    boot_results <- replicate(bootstrap_n, {
      sample_idx <- sample(1:nrow(boot_data), replace = TRUE)
      boot_fun(boot_data, sample_idx)
    }, simplify = FALSE)

    # Convert to matrix and handle failures
    boot_matrix <- do.call(cbind, boot_results)
    rownames(boot_matrix) <- model_coef_names_no_intercept

    # Check if we have any successful bootstrap samples
    successful_boots <- apply(boot_matrix, 2, function(x) !all(is.na(x)))
    n_successful_boot <- sum(successful_boots)

    if (verbose) cat(paste("Successful bootstrap samples:", n_successful_boot, "out of", bootstrap_n, "\n"))

    if (n_successful_boot < 10) {
      warning("Very few successful bootstrap samples (", n_successful_boot, "). Results may be unreliable.")
    }

    # Calculate bootstrap summary only from successful samples
    if (n_successful_boot > 0) {
      bootstrap_summary <- tibble(
        term = model_coef_names_no_intercept,
        estimate = rowMeans(boot_matrix[, successful_boots, drop = FALSE], na.rm = TRUE),
        se = apply(boot_matrix[, successful_boots, drop = FALSE], 1, sd, na.rm = TRUE),
        ci_lower = apply(boot_matrix[, successful_boots, drop = FALSE], 1, quantile, probs = 0.025, na.rm = TRUE),
        ci_upper = apply(boot_matrix[, successful_boots, drop = FALSE], 1, quantile, probs = 0.975, na.rm = TRUE)
      )

      bootstrap_results <- boot_matrix

      if (verbose) {
        cat("\nBootstrap Results (from imputation", bootstrap_imputation, "):\n")
        print(bootstrap_summary, digits = 4)
      }
    } else {
      # Fallback if no bootstrap samples succeeded
      warning("No successful bootstrap samples. Using GLM standard errors.")
      bootstrap_summary <- tibble(
        term = model_coef_names_no_intercept,
        estimate = coef(imputation_results[[bootstrap_imputation]]$model)[model_coef_names_no_intercept],
        se = summary(imputation_results[[bootstrap_imputation]]$model)$coefficients[model_coef_names_no_intercept, "Std. Error"],
        ci_lower = NA_real_,
        ci_upper = NA_real_
      )
      bootstrap_results <- boot_matrix
    }
  }

  # 7. BALANCE SUMMARY
  if (verbose) cat("\nStep 7: Summarizing balance across imputations...\n")

  balance_summary <- map_dfr(1:length(imputation_results), function(i) {
    bal <- imputation_results[[i]]$balance
    if ("Balance" %in% names(bal)) {
      balance_df <- as.data.frame(bal$Balance)

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

  final_n <- median(sapply(imputation_results, function(x) x$n))

  if (verbose) cat(paste("\nAnalysis completed. Median sample size across imputations:", final_n, "\n"))

  # 8. RETURN RESULTS (maintaining naming compatibility)
  return(list(
    # Original outputs (from first imputation for compatibility)
    data = imputation_results[[1]]$data,
    model = imputation_results[[1]]$model,
    weights = imputation_results[[1]]$weights,
    balance = imputation_results[[1]]$balance,
    bootstrap_summary = bootstrap_summary,  # Bootstrap from ONE imputation
    bootstrap_results = bootstrap_results,  # Bootstrap from ONE imputation
    sample_sizes = list(
      initial = initial_n,
      final = final_n
    ),

    # New MI-specific outputs
    pooled_results = pooled_results,        # PRIMARY RESULTS - use these!
    imputation_results = imputation_results,
    balance_summary = balance_summary,
    n_imputations = m,
    n_bootstraps = bootstrap_n,
    bootstrap_imputation_used = bootstrap_imputation
  ))
}
