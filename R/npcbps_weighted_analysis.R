#' Nonparametric CBPS Propensity Score Weighted GLM Analysis with Proper Multiple Imputation
#'
#' Performs comprehensive Nonparametric Covariate Balancing Propensity Score (NPCBPS) weighted analysis
#' with proper multiple imputation using Rubin's rules, plus bootstrap of a single outcome model.
#'
#' REVISION NOTES (changes from previous version):
#' \itemize{
#'   \item NEW: \code{family} now accepts any glm family (object, function, or name,
#'     e.g. \code{quasibinomial()}, \code{quasibinomial}, or \code{"quasibinomial"})
#'     or the string \code{"ordinal"}, which fits a proportional odds model with
#'     \code{MASS::polr()}.
#'   \item NEW: Average marginal effects (AMEs) of the treatment are computed in
#'     every imputation with \code{marginaleffects}, weighted by the NPCBPS weights,
#'     and pooled with Rubin's rules (\code{pooled_ame}). For ordinal outcomes the
#'     default AME is on the expected category score (scores 0, 1, 2, ...);
#'     \code{ame_hypothesis} lets you supply other weights, e.g. c(0, 0, 0, 1, 1)
#'     for the AME on P(top two categories).
#'   \item NEW: \code{ame_scale} multiplies the AME and its SE, e.g. 2 to report a
#'     fractional-logit AME for an outcome rescaled as (y - 1) / 2 in original units.
#'   \item FIXED: ordinal models pool only the slope coefficients (thresholds are
#'     excluded), so coefficient and variance vectors stay aligned.
#'   \item FIXED: convergence checks use \code{isTRUE()}, so models without a
#'     \code{$converged} element no longer fail silently in the bootstrap.
#'   \item FIXED: Barnard-Rubin degrees of freedom no longer return NaN when the
#'     between-imputation variance is zero (e.g. no missing data).
#'   \item MASS and marginaleffects are used via \code{::} and are not attached,
#'     so MASS::select() does not mask dplyr::select() in your session.
#' }
#'
#' Earlier revision notes:
#' \itemize{
#'   \item Balance metric depends on treatment type: |SMD| plus variance ratios for
#'     binary treatments, |r| for continuous treatments.
#'   \item Balance is collected from ALL imputations and summarized per covariate.
#'   \item Pooled standard errors use robust (sandwich) variance estimates by default.
#' }
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
#' @param cbps_estimand Character. CBPS estimand: "ATE" or "ATT" (default: "ATE")
#' @param bootstrap_n Integer. Number of bootstrap samples for the outcome model (default: 1000, set to 0 to skip)
#' @param bootstrap_seed Integer. Random seed for bootstrap (default: 20250417)
#' @param bootstrap_imputation Integer. Which imputation to use for bootstrap (default: 1)
#' @param balance_threshold_m Numeric. |SMD| threshold for binary treatments (default: 0.1)
#' @param balance_threshold_v Numeric. Variance ratio threshold for binary treatments (default: 2)
#' @param balance_threshold_r Numeric. |r| threshold for continuous treatments (default: 0.1)
#' @param robust_se Logical. Use robust (sandwich) standard errors for pooling (default: TRUE)
#' @param se_type Character. Sandwich estimator type passed to \code{sandwich::vcovHC} for glm
#'   families (default: "HC0"). Ordinal models always use \code{sandwich::sandwich()} (HC0).
#' @param family A glm family (object, function, or name) or the string "ordinal"
#'   (default: gaussian()). For binomial/quasibinomial the outcome must lie in [0, 1].
#'   For "ordinal" the outcome is converted to an ordered factor and must have at least 3 levels.
#' @param ame_hypothesis Numeric vector or NULL. Ordinal outcomes only: weights applied to the
#'   per-category AMEs, one per outcome level in order. NULL (default) uses scores 0, 1, 2, ...
#'   (AME on the expected category score).
#' @param ame_scale Numeric. Multiplier applied to the AME and its SE (default: 1).
#' @param verbose Logical. Whether to print progress messages (default: TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{data}, \code{model}, \code{weights}, \code{balance}: from the first imputation (compatibility)
#'   \item \code{balance_by_imputation}: balance statistics for every covariate in every imputation
#'   \item \code{balance_summary}: per-covariate balance across imputations (USE THIS FOR REPORTING)
#'   \item \code{treatment_type}: "binary" or "continuous"
#'   \item \code{family_label}: description of the outcome model
#'   \item \code{pooled_ame}: Rubin's rules pooled average marginal effect of the treatment (PRIMARY EFFECT)
#'   \item \code{ame_by_imputation}: AME estimate and variance from each imputation
#'   \item \code{pooled_results}: Rubin's rules pooled model coefficients (link scale for non-gaussian families)
#'   \item \code{bootstrap_summary}, \code{bootstrap_results}: bootstrap of the selected imputation's coefficients
#'   \item \code{imputation_results}: list of results from each imputation
#'   \item \code{sample_sizes}: initial and final sample sizes
#' }
#'
#' @details
#' Workflow:
#' \enumerate{
#'   \item Data cleaning and missing data assessment with MCAR testing
#'   \item Multiple imputation using MICE (m imputations kept separate)
#'   \item For EACH imputation: NPCBPS weights, balance assessment, weighted outcome
#'     model with robust SEs, and the weighted AME of the treatment
#'   \item Summarize balance across all imputations
#'   \item Pool coefficients and AMEs using Rubin's rules
#'   \item Bootstrap ONE of the fitted outcome models
#' }
#'
#' AMEs: for a continuous treatment, the weighted average derivative of the expected outcome
#' with respect to the treatment (\code{marginaleffects::avg_slopes}); for a binary treatment,
#' the weighted average contrast between its two values (\code{marginaleffects::avg_comparisons}).
#' Estimates are on the outcome's response scale. For a gaussian model with no interactions
#' the AME equals the treatment coefficient, which is a useful check.
#'
#' Robust standard errors treat the estimated weights as fixed. They correct the main problem
#' with model-based standard errors under propensity weighting, but do not account for
#' uncertainty in estimating the weights themselves.
#'
#' @references
#' Allison, P. (2015). Imputation by predictive mean matching: Promise & peril. Statistical Horizons.
#'
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
#' continuous treatment: Application to the efficacy of political advertisements.
#' The Annals of Applied Statistics, 12(1), 156-177.
#'
#' Papke, L. E., & Wooldridge, J. M. (1996). Econometric methods for fractional response
#' variables with an application to 401(k) plan participation rates. Journal of Applied
#' Econometrics, 11(6), 619-632.
#'
#' Rubin, D. B. (1987). Multiple Imputation for Nonresponse in Surveys. Wiley.
#'
#' van Buuren, S. (2018). Flexible Imputation of Missing Data (2nd ed.). CRC Press.
#'
#' @examples
#' \dontrun{
#' # Fractional logit for mpss rescaled to [0, 1]; AME reported in mpss units
#' res_mpss <- npcbps_weighted_analysis(
#'   data = mad_edits, outcome_var = "mpss.mgm_rescale", treatment_var = "my_treatment",
#'   imputation_vars = c("mother_age"), family = quasibinomial(), ame_scale = 2
#' )
#'
#' # Ordinal pregnancy advice (expected-score AME)
#' res_adv <- npcbps_weighted_analysis(
#'   data = mad_edits, outcome_var = "preg_ad_ord", treatment_var = "my_treatment",
#'   imputation_vars = c("mother_age"), family = "ordinal"
#' )
#'
#' # Ordinal food provisioning, AME on P(more than once a week or every day)
#' res_food <- npcbps_weighted_analysis(
#'   data = mad_edits, outcome_var = "food_ord", treatment_var = "my_treatment",
#'   imputation_vars = c("mother_age"), family = "ordinal",
#'   ame_hypothesis = c(0, 0, 0, 1, 1)
#' )
#'
#' res_food$pooled_ame        # PRIMARY effect estimate
#' res_food$pooled_results    # model coefficients (log-odds)
#' res_food$balance_summary   # balance across imputations
#' }
#'
#' @export
#' @importFrom dplyr select all_of filter if_all mutate bind_cols bind_rows summarise across everything group_by case_when if_else n
#' @importFrom mice mice complete
#' @importFrom purrr map map_dfr
#' @importFrom WeightIt weightit
#' @importFrom cobalt bal.tab
#' @importFrom tibble tibble
#' @importFrom sandwich vcovHC sandwich
#' @importFrom stats glm as.formula coef quantile sd vcov qt pt var setNames
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
    balance_threshold_r = 0.1,
    robust_se = TRUE,
    se_type = "HC0",
    family = gaussian(),
    ame_hypothesis = NULL,
    ame_scale = 1,
    verbose = TRUE
) {

  # ---- NEW: normalize the family argument ----
  is_ordinal <- is.character(family) && length(family) == 1 && tolower(family) == "ordinal"
  if (!is_ordinal) {
    if (is.character(family)) family <- get(family, mode = "function")
    if (is.function(family))  family <- family()
    if (!inherits(family, "family")) {
      stop("`family` must be a glm family (e.g. gaussian(), quasibinomial()) or \"ordinal\".")
    }
  }
  fam_label <- if (is_ordinal) {
    "ordinal (proportional odds, logit link)"
  } else {
    paste(family$family, "with", family$link, "link")
  }
  if (!is.null(ame_hypothesis) && !is_ordinal) {
    warning("`ame_hypothesis` only applies to ordinal outcomes and will be ignored.", call. = FALSE)
  }

  # Load required libraries (attached)
  required_packages <- c("dplyr", "mice", "purrr", "WeightIt", "cobalt", "tibble")
  if (robust_se) required_packages <- c(required_packages, "sandwich")
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  }
  # ---- NEW: packages used via :: only (not attached; avoids MASS::select masking dplyr) ----
  namespace_packages <- c("marginaleffects", if (is_ordinal) "MASS")
  for (pkg in namespace_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
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

  # ---- NEW: prepare and check the outcome for the chosen family ----
  if (is_ordinal) {
    df_clean[[outcome_var]] <- droplevels(factor(df_clean[[outcome_var]], ordered = TRUE))
    n_out_levels <- nlevels(df_clean[[outcome_var]])
    if (n_out_levels < 3) {
      stop("Ordinal outcome '", outcome_var, "' has only ", n_out_levels,
           " observed levels; polr needs at least 3. Use family = quasibinomial() for a binary outcome.")
    }
    if (!is.null(ame_hypothesis) && length(ame_hypothesis) != n_out_levels) {
      stop("`ame_hypothesis` has length ", length(ame_hypothesis), " but '", outcome_var,
           "' has ", n_out_levels, " observed levels (",
           paste(levels(df_clean[[outcome_var]]), collapse = ", "), ").")
    }
    if (verbose) {
      cat("Ordinal outcome levels:", paste(levels(df_clean[[outcome_var]]), collapse = " < "), "\n")
      cat("Level counts:", paste(table(df_clean[[outcome_var]]), collapse = ", "), "\n")
    }
  } else if (family$family %in% c("binomial", "quasibinomial")) {
    y <- df_clean[[outcome_var]]
    if (!is.numeric(y) || any(y < 0 | y > 1)) {
      stop("Outcome '", outcome_var, "' must be numeric in [0, 1] for ", family$family,
           " (range observed: ", paste(round(range(y), 3), collapse = " to "),
           "). Rescale it first, e.g. (y - min) / (max - min).")
    }
  }

  # Determine treatment type (drives which balance metric is used)
  treat_vals <- df_clean[[treatment_var]]
  n_treat_levels <- length(unique(treat_vals))

  if (n_treat_levels == 2) {
    treatment_type <- "binary"
  } else if (is.numeric(treat_vals)) {
    treatment_type <- "continuous"
  } else {
    stop("Treatment '", treatment_var, "' has more than two categories. ",
         "Balance summaries currently support binary or continuous treatments only.")
  }

  balance_metric <- if (treatment_type == "continuous") "r" else "SMD"
  threshold_primary <- if (treatment_type == "continuous") balance_threshold_r else balance_threshold_m
  vr_low <- 1 / balance_threshold_v
  vr_high <- balance_threshold_v

  if (verbose) {
    cat("Treatment type detected:", treatment_type, "\n")
    cat("Primary balance metric: |", balance_metric, "| < ", threshold_primary, "\n", sep = "")
    if (treatment_type == "binary") {
      cat("Secondary balance metric: VR between", vr_low, "and", vr_high, "\n")
    }
  }

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
            cat("     No evidence against MCAR (Little's test not significant)\n")
            cat("     Multiple imputation assumptions are reasonable\n")
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

  # Helper to pull the correct balance statistics from a cobalt table
  extract_balance <- function(bal, imp) {
    b <- as.data.frame(bal$Balance)

    needed <- if (treatment_type == "continuous") c("Corr.Un", "Corr.Adj") else c("Diff.Un", "Diff.Adj")
    missing_cols <- setdiff(needed, names(b))
    if (length(missing_cols) > 0) {
      stop("Expected balance columns not found in bal.tab output: ",
           paste(missing_cols, collapse = ", "),
           ". Available columns: ", paste(names(b), collapse = ", "))
    }

    tibble::tibble(
      imputation     = imp,
      covariate      = rownames(b),
      covariate_type = as.character(b$Type),
      metric         = balance_metric,
      before         = b[[needed[1]]],
      after          = b[[needed[2]]],
      vr_before      = if ("V.Ratio.Un"  %in% names(b)) b[["V.Ratio.Un"]]  else NA_real_,
      vr_after       = if ("V.Ratio.Adj" %in% names(b)) b[["V.Ratio.Adj"]] else NA_real_
    )
  }

  bal_stats <- if (treatment_type == "continuous") {
    "correlations"
  } else {
    c("mean.diffs", "variance.ratios")
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
    cat("DEBUG: Outcome family:", fam_label, "\n")
    cat("DEBUG: Standard errors:",
        if (!robust_se) "model-based" else if (is_ordinal) "robust (sandwich, HC0)" else paste0("robust (", se_type, ")"),
        "\n")
  }

  # ---- NEW: single fitting function for glm families and ordinal ----
  fit_outcome <- function(d, w) {
    d$.w <- w
    if (is_ordinal) {
      d[[outcome_var]] <- droplevels(d[[outcome_var]])
      m <- MASS::polr(outcome_formula, data = d, weights = .w,
                      Hess = TRUE, method = "logistic")
      m$converged <- isTRUE(m$convergence == 0)
    } else {
      m <- glm(outcome_formula, data = d, weights = .w, family = family)
    }
    m
  }

  # ---- NEW: variance-covariance for pooling (full = TRUE keeps polr thresholds) ----
  vcov_fallback_warned <- FALSE
  get_vcov <- function(m, full = FALSE) {
    V <- if (!robust_se) {
      vcov(m)
    } else if (is_ordinal) {
      tryCatch(sandwich::sandwich(m), error = function(e) {
        if (!vcov_fallback_warned) {
          warning("Robust vcov failed for polr (", e$message,
                  "); using model-based vcov instead.", call. = FALSE)
          vcov_fallback_warned <<- TRUE
        }
        vcov(m)
      })
    } else {
      sandwich::vcovHC(m, type = se_type)
    }
    if (full) return(V)
    keep <- names(coef(m))
    V[keep, keep, drop = FALSE]
  }

  # ---- NEW: weighted average marginal effect of the treatment ----
  compute_ame <- function(m, d, w) {
    d$.w <- w
    args <- list(m, newdata = d, wts = w, vcov = get_vcov(m, full = TRUE))

    if (treatment_type == "binary") {
      tv <- d[[treatment_var]]
      contrast_vals <- if (is.factor(tv)) levels(droplevels(tv)) else sort(unique(tv))
      args$variables <- stats::setNames(list(contrast_vals), treatment_var)
      me_fun <- marginaleffects::avg_comparisons
    } else {
      args$variables <- treatment_var
      me_fun <- marginaleffects::avg_slopes
    }

    ame <- tryCatch(
      do.call(me_fun, args),
      error = function(e) {
        # Retry with model-based vcov if the robust matrix is not accepted
        args$vcov <- TRUE
        warning("AME with robust vcov failed (", e$message,
                "); retried with model-based vcov.", call. = FALSE)
        do.call(me_fun, args)
      }
    )

    if (is_ordinal) {
      lv <- levels(droplevels(d[[outcome_var]]))
      ord <- match(lv, as.character(ame$group))
      if (anyNA(ord)) stop("Could not match AME rows to outcome levels.")
      est <- ame$estimate[ord]
      V_ame <- as.matrix(stats::vcov(ame))[ord, ord, drop = FALSE]
      h <- if (is.null(ame_hypothesis)) seq_along(lv) - 1 else ame_hypothesis
      list(estimate = sum(h * est) * ame_scale,
           variance = drop(t(h) %*% V_ame %*% h) * ame_scale^2)
    } else {
      if (nrow(ame) != 1) stop("Expected one AME row, got ", nrow(ame), ".")
      list(estimate = ame$estimate * ame_scale,
           variance = ame$std.error^2 * ame_scale^2)
    }
  }

  imputation_results <- list()

  for (imp in 1:mice_m) {
    if (verbose) cat(paste("  Processing imputation", imp, "of", mice_m, "...\n"))

    df_imp <- bind_cols(outcome_treatment_data, imputed_datasets[[imp]])

    missing_check <- df_imp %>%
      dplyr::select(all_of(ps_covariate_vars)) %>%
      summarise(across(everything(), ~sum(is.na(.x))))

    if (any(missing_check > 0)) {
      df_imp <- df_imp %>%
        filter(if_all(all_of(ps_covariate_vars), ~!is.na(.x)))
      if (verbose) cat(paste("    Removed", initial_n - nrow(df_imp), "rows with missing covariates\n"))
    }

    if (verbose) cat("    Note: NPCBPS uses kernel-based estimation - may take longer than parametric CBPS\n")

    tryCatch({
      npcbps_weights <- WeightIt::weightit(
        ps_formula,
        data = df_imp,
        method = "npcbps",
        estimand = cbps_estimand
      )

      balance_table <- cobalt::bal.tab(npcbps_weights,
                                       stats = bal_stats,
                                       un = TRUE)

      balance_long <- extract_balance(balance_table, imp)

      if (verbose) {
        cat(paste0("    Balance: max |", balance_metric, "| after weighting = ",
                   round(max(abs(balance_long$after), na.rm = TRUE), 4), "\n"))
      }

      # ---- CHANGED: family-aware fit and vcov ----
      weighted_model <- fit_outcome(df_imp, npcbps_weights$weights)
      model_vcov <- get_vcov(weighted_model)

      # ---- NEW: AME (failure here does not drop the imputation) ----
      ame_i <- tryCatch(
        compute_ame(weighted_model, df_imp, npcbps_weights$weights),
        error = function(e) {
          warning(paste("AME failed in imputation", imp, ":", e$message), call. = FALSE)
          list(estimate = NA_real_, variance = NA_real_)
        }
      )

      imputation_results[[imp]] <- list(
        data = df_imp,
        weights = npcbps_weights,
        balance = balance_table,
        balance_long = balance_long,
        model = weighted_model,
        coefficients = coef(weighted_model),
        vcov = model_vcov,
        ame_estimate = ame_i$estimate,
        ame_variance = ame_i$variance,
        n = nrow(df_imp)
      )

      if (verbose) {
        cat(paste("    n =", nrow(df_imp), "- Model converged:", isTRUE(weighted_model$converged),
                  "- AME =", round(ame_i$estimate, 4), "\n"))
      }

    }, error = function(e) {
      warning(paste("Error in imputation", imp, ":", e$message), call. = FALSE)
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

  # 5. BALANCE ACROSS ALL IMPUTATIONS
  if (verbose) cat("Step 5: Summarizing balance across all imputations...\n")

  balance_by_imputation <- dplyr::bind_rows(lapply(imputation_results, `[[`, "balance_long"))

  safe_min  <- function(x) if (all(is.na(x))) NA_real_ else min(x, na.rm = TRUE)
  safe_max  <- function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE)
  safe_mean <- function(x) if (all(is.na(x))) NA_real_ else mean(x, na.rm = TRUE)

  balance_summary <- balance_by_imputation %>%
    dplyr::group_by(covariate, covariate_type, metric) %>%
    dplyr::summarise(
      mean_abs_before = safe_mean(abs(before)),
      mean_abs_after  = safe_mean(abs(after)),
      max_abs_after   = safe_max(abs(after)),
      mean_vr_after   = safe_mean(vr_after),
      min_vr_after    = safe_min(vr_after),
      max_vr_after    = safe_max(vr_after),
      n_imputations   = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      primary_balanced = max_abs_after < threshold_primary,
      vr_balanced = dplyr::if_else(
        is.na(min_vr_after), NA,
        min_vr_after >= vr_low & max_vr_after <= vr_high
      ),
      balance_status = dplyr::case_when(
        !primary_balanced   ~ paste0("Imbalanced: |", metric, "| >= ", threshold_primary),
        vr_balanced %in% FALSE ~ "Primary met; VR outside range",
        TRUE                ~ paste0("Balanced: |", metric, "| < ", threshold_primary)
      )
    )

  primary_fail <- balance_summary$covariate[!balance_summary$primary_balanced]
  vr_fail      <- balance_summary$covariate[balance_summary$vr_balanced %in% FALSE]

  if (length(primary_fail) > 0) {
    warning(paste0("Poor balance (|", balance_metric, "| >= ", threshold_primary,
                   " in at least one imputation) for: ",
                   paste(primary_fail, collapse = ", ")),
            call. = FALSE)
  }
  if (length(vr_fail) > 0) {
    warning(paste0("Variance ratio outside [", vr_low, ", ", vr_high,
                   "] in at least one imputation for: ",
                   paste(vr_fail, collapse = ", ")),
            call. = FALSE)
  }

  if (verbose) {
    cat("  Balance summary across", n_successful, "imputations:\n")
    print(balance_summary %>%
            dplyr::select(covariate, metric, mean_abs_after, max_abs_after,
                          min_vr_after, max_vr_after, balance_status),
          digits = 4)
  }

  # Get coefficient names (polr has no intercept; its thresholds are not in coef())
  all_coef_names <- names(imputation_results[[1]]$coefficients)
  model_coef_names_no_intercept <- all_coef_names[all_coef_names != "(Intercept)"]

  if (verbose) {
    cat("DEBUG: Model coefficient names:", paste(all_coef_names, collapse = ", "), "\n")
    cat("DEBUG: Coefficients for bootstrap:", paste(model_coef_names_no_intercept, collapse = ", "), "\n")
  }

  # 6. POOL RESULTS USING RUBIN'S RULES
  if (verbose) cat("Step 6: Pooling results using Rubin's rules...\n")

  # ---- NEW: reusable Rubin's rules (rows = parameters, columns = imputations) ----
  pool_rubin <- function(Q_m, U_m, df_obs) {
    m <- ncol(Q_m)
    Q_bar <- rowMeans(Q_m)
    U_bar <- rowMeans(U_m)
    B <- if (m > 1) apply(Q_m, 1, stats::var) else rep(0, nrow(Q_m))
    T_var <- U_bar + B + B / m
    SE <- sqrt(T_var)

    # Degrees of freedom (Barnard & Rubin, 1999); falls back to df_adj when B = 0
    lambda <- (B + B / m) / T_var
    df_old <- (m - 1) / lambda^2
    df_adj <- (df_obs + 1) / (df_obs + 3) * df_obs * (1 - lambda)
    df <- ifelse(is.finite(df_old), (df_old * df_adj) / (df_old + df_adj), df_adj)

    t_stat <- Q_bar / SE
    r <- (B + B / m) / U_bar
    tibble(
      term = rownames(Q_m),
      estimate = Q_bar,
      se = SE,
      statistic = t_stat,
      p.value = 2 * pt(abs(t_stat), df, lower.tail = FALSE),
      ci_lower = Q_bar - qt(0.975, df) * SE,
      ci_upper = Q_bar + qt(0.975, df) * SE,
      df = df,
      fmi = (r + 2 / (df + 3)) / (r + 1),
      within_var = U_bar,
      between_var = B,
      total_var = T_var
    )
  }

  m <- n_successful
  df_obs <- imputation_results[[1]]$model$df.residual
  se_label <- if (!robust_se) "model_based" else if (is_ordinal) "robust_sandwich_HC0" else paste0("robust_", se_type)

  # Coefficients (cbind keeps a parameters x imputations matrix even with one coefficient)
  Q_m <- do.call(cbind, lapply(imputation_results, function(x) x$coefficients))
  U_m <- do.call(cbind, lapply(imputation_results, function(x) diag(x$vcov)))
  rownames(Q_m) <- rownames(U_m) <- all_coef_names

  pooled_results <- pool_rubin(Q_m, U_m, df_obs) %>%
    dplyr::mutate(se_method = se_label)

  # ---- NEW: pooled AME ----
  ame_by_imputation <- tibble(
    imputation = seq_len(m),
    estimate = sapply(imputation_results, `[[`, "ame_estimate"),
    variance = sapply(imputation_results, `[[`, "ame_variance")
  )
  ame_ok <- !is.na(ame_by_imputation$estimate) & !is.na(ame_by_imputation$variance)

  pooled_ame <- NULL
  if (sum(ame_ok) > 0) {
    if (sum(ame_ok) < m) {
      warning(paste("AME available for only", sum(ame_ok), "of", m, "imputations; pooling those."),
              call. = FALSE)
    }
    ame_label <- paste0("AME_", treatment_var)
    Q_ame <- matrix(ame_by_imputation$estimate[ame_ok], nrow = 1, dimnames = list(ame_label, NULL))
    U_ame <- matrix(ame_by_imputation$variance[ame_ok], nrow = 1, dimnames = list(ame_label, NULL))
    pooled_ame <- pool_rubin(Q_ame, U_ame, df_obs) %>%
      dplyr::mutate(
        se_method = se_label,
        ame_type = if (is_ordinal) {
          if (is.null(ame_hypothesis)) "expected category score" else
            paste0("weighted categories (", paste(ame_hypothesis, collapse = ", "), ")")
        } else "response scale",
        ame_scale = ame_scale,
        n_imputations = sum(ame_ok)
      )
  } else {
    warning("AME could not be computed in any imputation; see earlier warnings.", call. = FALSE)
  }

  if (verbose) {
    cat("\nPooled Results (Rubin's Rules - PRIMARY INFERENCE):\n")
    cat("Number of imputations:", m, "\n")
    cat("Outcome model:", fam_label, "\n")
    if (!is.null(pooled_ame)) {
      cat("\nAverage marginal effect of", treatment_var, "(", pooled_ame$ame_type, "):\n")
      print(pooled_ame %>%
              dplyr::select(term, estimate, se, ci_lower, ci_upper, p.value, fmi),
            digits = 4)
    }
    cat("\nCoefficients", if (is_ordinal || family$family != "gaussian") "(link scale)" else "", ":\n")
    print(pooled_results %>%
            dplyr::select(term, estimate, se, ci_lower, ci_upper, p.value, fmi) %>%
            dplyr::filter(term != "(Intercept)"),
          digits = 4)
  }

  # 7. BOOTSTRAP ONE OUTCOME MODEL
  bootstrap_summary <- NULL
  bootstrap_results <- NULL

  if (bootstrap_n > 0) {
    if (verbose) cat("\nStep 7: Bootstrapping outcome model from imputation", bootstrap_imputation, "...\n")

    boot_data <- imputation_results[[bootstrap_imputation]]$data %>%
      mutate(weight = imputation_results[[bootstrap_imputation]]$weights$weights)

    set.seed(bootstrap_seed)

    boot_fun <- function(data, indices) {
      tryCatch({
        d <- data[indices, ]

        if (length(unique(d[[treatment_var]])) < 2) {
          return(rep(NA, length(model_coef_names_no_intercept)))
        }

        # ---- CHANGED: family-aware fit; isTRUE() so missing $converged can't error ----
        model <- fit_outcome(d, d$weight)

        if (!isTRUE(model$converged)) {
          return(rep(NA, length(model_coef_names_no_intercept)))
        }

        result <- coef(model)[model_coef_names_no_intercept]

        if (length(result) != length(model_coef_names_no_intercept) || anyNA(result)) {
          return(rep(NA, length(model_coef_names_no_intercept)))
        }

        return(as.numeric(result))

      }, error = function(e) {
        return(rep(NA, length(model_coef_names_no_intercept)))
      })
    }

    boot_results <- replicate(bootstrap_n, {
      sample_idx <- sample(1:nrow(boot_data), replace = TRUE)
      boot_fun(boot_data, sample_idx)
    }, simplify = FALSE)

    boot_matrix <- do.call(cbind, boot_results)
    if (is.null(dim(boot_matrix))) boot_matrix <- matrix(boot_matrix, nrow = 1)
    rownames(boot_matrix) <- model_coef_names_no_intercept

    successful_boots <- apply(boot_matrix, 2, function(x) !all(is.na(x)))
    n_successful_boot <- sum(successful_boots)

    if (verbose) cat(paste("Successful bootstrap samples:", n_successful_boot, "out of", bootstrap_n, "\n"))

    if (n_successful_boot < 10) {
      warning("Very few successful bootstrap samples (", n_successful_boot, "). Results may be unreliable.")
    }

    if (n_successful_boot > 0) {
      bootstrap_summary <- tibble(
        term = model_coef_names_no_intercept,
        estimate = rowMeans(boot_matrix[, successful_boots, drop = FALSE], na.rm = TRUE),
        se = apply(boot_matrix[, successful_boots, drop = FALSE], 1, sd, na.rm = TRUE),
        ci_lower = apply(boot_matrix[, successful_boots, drop = FALSE], 1, quantile, probs = 0.025, na.rm = TRUE),
        ci_upper = apply(boot_matrix[, successful_boots, drop = FALSE], 1, quantile, probs = 0.975, na.rm = TRUE),
        n_successful = n_successful_boot
      )

      bootstrap_results <- boot_matrix

      if (verbose) {
        cat("\nBootstrap Results (from imputation", bootstrap_imputation, "):\n")
        print(bootstrap_summary, digits = 4)
      }
    } else {
      warning("No successful bootstrap samples. Using model standard errors.")
      bootstrap_summary <- tibble(
        term = model_coef_names_no_intercept,
        estimate = coef(imputation_results[[bootstrap_imputation]]$model)[model_coef_names_no_intercept],
        se = sqrt(diag(imputation_results[[bootstrap_imputation]]$vcov))[model_coef_names_no_intercept],
        ci_lower = NA_real_,
        ci_upper = NA_real_,
        n_successful = 0L
      )
      bootstrap_results <- boot_matrix
    }
  }

  final_n <- median(sapply(imputation_results, function(x) x$n))

  if (verbose) cat(paste("\nAnalysis completed. Median sample size across imputations:", final_n, "\n"))

  # 8. RETURN RESULTS
  return(list(
    data = imputation_results[[1]]$data,
    model = imputation_results[[1]]$model,
    weights = imputation_results[[1]]$weights,
    balance = imputation_results[[1]]$balance,
    bootstrap_summary = bootstrap_summary,
    bootstrap_results = bootstrap_results,
    sample_sizes = list(
      initial = initial_n,
      final = final_n
    ),

    pooled_ame = pooled_ame,                       # PRIMARY EFFECT ESTIMATE (outcome scale)
    ame_by_imputation = ame_by_imputation,
    pooled_results = pooled_results,               # coefficients (link scale if non-gaussian)
    imputation_results = imputation_results,
    balance_by_imputation = balance_by_imputation,
    balance_summary = balance_summary,
    treatment_type = treatment_type,
    family_label = fam_label,
    n_imputations = m,
    n_bootstraps = bootstrap_n,
    bootstrap_imputation_used = bootstrap_imputation
  ))
}
