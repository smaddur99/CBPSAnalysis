#' Check Treatment Significance Across Multiple Models
#'
#' Analyzes treatment significance across multiple CBPS model specifications,
#' focusing specifically on the primary treatment variable of interest.
#'
#' @param model_results_list List of CBPS analysis results from \code{\link{cbps_weighted_analysis}}, or a single result object
#' @param model_names Character vector of names for the models (optional)
#' @param primary_treatment Character. Name of the primary treatment variable to focus on
#' @param additional_predictors Character vector. Additional predictor names (used for filtering during processing)
#' @param alpha_levels Numeric vector. Alpha levels for significance testing (default: c(0.05, 0.01, 0.001))
#' @param include_bootstrap Logical. Whether to include bootstrap significance tests (default: TRUE)
#' @param exclude_intercept Logical. Whether to exclude intercept terms from results (default: TRUE)
#' @param focus_on_treatment_only Logical. Whether to show only the primary treatment variable (default: TRUE)
#' @param verbose Logical. Whether to print detailed summary information (default: TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{detailed_results}: Full results for all variables and models
#'   \item \code{summary_table}: Summary statistics for each model
#'   \item \code{significant_only}: Results filtered to significant findings only
#'   \item \code{alpha_levels}: Alpha levels used for testing
#'   \item \code{n_models}: Number of models analyzed
#'   \item \code{model_names}: Names of the models
#' }
#'
#' @details
#' This function is designed to help researchers compare treatment effects across
#' different model specifications. It provides:
#' \itemize{
#'   \item Significance testing at multiple alpha levels
#'   \item Bootstrap confidence interval assessment
#'   \item Balance quality metrics for each model
#'   \item Model fit statistics (AIC, BIC, BICc)
#'   \item Publication-ready summaries
#' }
#'
#' When \code{focus_on_treatment_only = TRUE}, the function filters results to show
#' only the primary treatment variable, making it easier to compare treatment effects
#' across different model specifications.
#'
#' The function automatically handles balance assessment, providing warnings when
#' covariate balance is poor and might affect the reliability of results.
#'
#' @examples
#' \dontrun{
#' # Run multiple CBPS models
#' model1 <- cbps_weighted_analysis(data, "outcome", "treatment",
#'                                  imputation_vars = c("age", "income"))
#' model2 <- cbps_weighted_analysis(data, "outcome", "treatment",
#'                                  imputation_vars = c("age", "income", "education"))
#'
#' # Check significance across models
#' sig_results <- check_treatment_significance(
#'   model_results_list = list(model1, model2),
#'   model_names = c("Basic Model", "Extended Model"),
#'   primary_treatment = "treatment"
#' )
#'
#' # View summary
#' sig_results$summary_table
#'
#' # View significant results only
#' sig_results$significant_only
#' }
#'
#' @seealso \code{\link{cbps_weighted_analysis}}, \code{\link{quick_sig_check}}
#' @export
#' @importFrom dplyr bind_rows group_by summarise across starts_with filter arrange slice_min select
#' @importFrom tibble tibble
#' @importFrom broom tidy
check_treatment_significance <- function(
    model_results_list,
    model_names = NULL,
    primary_treatment = "treatment_var",
    additional_predictors = NULL,
    alpha_levels = c(0.05, 0.01, 0.001),
    include_bootstrap = TRUE,
    exclude_intercept = TRUE,
    focus_on_treatment_only = TRUE,
    use_pooled = TRUE,
    verbose = TRUE
) {

  # Load required libraries
  required_packages <- c("dplyr", "tibble")
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  }

  # Handle single model input (convert to list)
  if (!is.list(model_results_list) || "model" %in% names(model_results_list)) {
    model_results_list <- list(model_results_list)
  }

  # Create model names if not provided
  if (is.null(model_names)) {
    model_names <- paste0("Model_", 1:length(model_results_list))
  }

  # Initialize results storage
  all_results <- list()

  # Process each model
  for (i in 1:length(model_results_list)) {

    model_result <- model_results_list[[i]]
    model_name <- model_names[i]

    if (verbose) cat("Processing", model_name, "...\n")

    # Initialize variables first (important for scope!)
    model_aic <- NA_real_
    model_bic <- NA_real_
    model_aicc <- NA_real_

    # Extract results using manual extraction
    tryCatch({

      # Check if input is from extract_cbps_results() or raw CBPS output
      if ("estimates" %in% names(model_result)) {
        # Input from extract_cbps_results() with include_balance_details = TRUE
        if (verbose) cat("  Using extracted CBPS results format\n")

        # Start with the estimates dataframe
        model_data <- model_result$estimates

        # Extract sample size from sample_sizes tibble if available
        if ("sample_sizes" %in% names(model_result)) {
          final_n <- model_result$sample_sizes %>%
            filter(metric == "final_n") %>%
            pull(value)
          model_data$sample_size <- final_n
        }

        # Extract balance statistics from summary_stats if available
        if ("summary_stats" %in% names(model_result)) {
          max_std_diff <- model_result$summary_stats %>%
            filter(metric == "max_abs_std_diff_after") %>%
            pull(value)
          mean_std_diff <- model_result$summary_stats %>%
            filter(metric == "mean_absolute_std_diff_after") %>%
            pull(value)
          n_balanced <- model_result$summary_stats %>%
            filter(metric == "n_variables_balanced") %>%
            pull(value)

          # Get total number of variables from balance data
          if ("balance" %in% names(model_result)) {
            n_total <- nrow(model_result$balance)
          } else {
            n_total <- NA_real_
          }

          model_data$max_std_diff <- max_std_diff
          model_data$mean_std_diff <- mean_std_diff
          model_data$n_balanced_vars <- n_balanced
          model_data$n_total_vars <- n_total
        } else {
          model_data$max_std_diff <- NA_real_
          model_data$mean_std_diff <- NA_real_
          model_data$n_balanced_vars <- NA_real_
          model_data$n_total_vars <- NA_real_
        }

        # Check for pooled results (MI analysis)
        if (!is.null(model_result$pooled_estimates) && use_pooled) {
          if (verbose) cat("  Using Rubin's rules pooled results\n")

          # Replace GLM estimates with pooled estimates where available
          pooled_data <- model_result$pooled_estimates

          model_data <- model_data %>%
            left_join(
              pooled_data %>%
                select(variable, pooled_estimate, pooled_se, pooled_statistic,
                       pooled_p_value, pooled_conf_low, pooled_conf_high,
                       pooled_df, pooled_fmi),
              by = "variable"
            ) %>%
            mutate(
              # Use pooled results for inference when available
              glm_estimate = ifelse(!is.na(pooled_estimate), pooled_estimate, glm_estimate),
              glm_se = ifelse(!is.na(pooled_se), pooled_se, glm_se),
              glm_statistic = ifelse(!is.na(pooled_statistic), pooled_statistic, glm_statistic),
              glm_p_value = ifelse(!is.na(pooled_p_value), pooled_p_value, glm_p_value),
              glm_conf_low = ifelse(!is.na(pooled_conf_low), pooled_conf_low, glm_conf_low),
              glm_conf_high = ifelse(!is.na(pooled_conf_high), pooled_conf_high, glm_conf_high),
              fmi = pooled_fmi,
              df_adjusted = pooled_df,
              inference_method = ifelse(!is.na(pooled_estimate), "Rubin's Rules (MI)", "Single Imputation")
            ) %>%
            select(-c(pooled_estimate, pooled_se, pooled_statistic, pooled_p_value,
                      pooled_conf_low, pooled_conf_high, pooled_df, pooled_fmi))
        } else {
          model_data <- model_data %>%
            mutate(
              inference_method = "Single Imputation",
              fmi = NA_real_,
              df_adjusted = NA_real_
            )
        }

        # Add model information
        model_data <- model_data %>%
          mutate(
            model_name = model_name,
            AIC = NA_real_,  # Not available from extracted results
            BIC = NA_real_,
            BICc = NA_real_,
            n_imputations = if(!is.null(model_result$pooled_estimates)) {
              # Infer from pooled results presence
              5  # Default assumption
            } else {
              1
            }
          )

      } else {
        # Original format: raw CBPS results
        if (verbose) cat("  Using raw CBPS results format\n")

        # Check if pooled results exist and should be used
        use_pooled_for_model <- use_pooled && !is.null(model_result$pooled_results)

        if (use_pooled_for_model) {
          # USE POOLED RESULTS (Rubin's rules)
          if (verbose) cat("  Using Rubin's rules pooled results\n")

          glm_summary <- model_result$pooled_results %>%
            rename(
              term = term,
              estimate = estimate,
              std.error = se,
              statistic = statistic,
              p.value = p.value,
              conf.low = ci_lower,
              conf.high = ci_upper
            ) %>%
            mutate(
              inference_method = "Rubin's Rules (MI)",
              fmi = fmi,
              df_adjusted = df
            )

        } else {
          # USE SINGLE IMPUTATION (backward compatible)
          if (verbose && !is.null(model_result$pooled_results)) {
            cat("  Using single imputation results (use_pooled=FALSE)\n")
          }

          glm_summary <- broom::tidy(model_result$model, conf.int = TRUE) %>%
            mutate(
              inference_method = "Single Imputation",
              fmi = NA_real_,
              df_adjusted = NA_real_
            )
        }

        # Check if confidence intervals were created
        if (!"conf.low" %in% names(glm_summary)) {
          glm_summary$conf.low <- NA_real_
          glm_summary$conf.high <- NA_real_
          if (verbose) cat("  Warning: No confidence intervals available for", model_name, "\n")
        }

        # Extract model fit statistics (AIC, BIC, BICc)
        tryCatch({
          model_aic <- AIC(model_result$model)
          model_bic <- BIC(model_result$model)

          # Calculate BICc (corrected BIC for small samples)
          n_obs <- nobs(model_result$model)
          k_params <- length(coef(model_result$model))

          model_aicc <- model_aic + (2 * k_params * (k_params + 1)) / (n_obs - k_params - 1)

          if (verbose) {
            cat("  Model fit statistics: AIC =", round(model_aic, 2),
                ", BIC =", round(model_bic, 2),
                ", BICc =", round(model_aicc, 2), "\n")
          }

        }, error = function(e) {
          if (verbose) cat("  Warning: Could not calculate AIC/BIC for", model_name, ":", e$message, "\n")
        })

        # Get bootstrap results
        bootstrap_summary <- model_result$bootstrap_summary

        # Combine results
        model_data <- glm_summary %>%
          dplyr::rename(
            variable = term,
            glm_estimate = estimate,
            glm_se = std.error,
            glm_statistic = statistic,
            glm_p_value = p.value,
            glm_conf_low = conf.low,
            glm_conf_high = conf.high
          ) %>%
          dplyr::left_join(
            bootstrap_summary %>%
              dplyr::rename(
                variable = term,
                bootstrap_estimate = estimate,
                bootstrap_se = se,
                bootstrap_conf_low = ci_lower,
                bootstrap_conf_high = ci_upper
              ),
            by = "variable"
          ) %>%
          dplyr::mutate(
            model_name = model_name,
            sample_size = model_result$sample_sizes$final,
            analysis_type = "Combined",
            AIC = model_aic,
            BIC = model_bic,
            BICc = model_aicc,
            n_imputations = if (!is.null(model_result$n_imputations)) model_result$n_imputations else 1,
            inference_method = first(inference_method)
          )

        # Initialize balance columns
        model_data$max_std_diff <- NA_real_
        model_data$mean_std_diff <- NA_real_
        model_data$n_balanced_vars <- NA_real_
        model_data$n_total_vars <- NA_real_

        # Extract overall balance statistics for the model
        if (!is.null(model_result$balance_summary) && use_pooled_for_model) {
          # Use average balance across imputations
          balance_sum <- model_result$balance_summary

          if (nrow(balance_sum) > 0) {
            max_std_diff <- mean(balance_sum$max_abs_std_diff, na.rm = TRUE)
            mean_std_diff <- mean(balance_sum$mean_abs_std_diff, na.rm = TRUE)
            n_balanced <- sum(balance_sum$n_above_threshold == 0, na.rm = TRUE)
            n_total <- nrow(balance_sum)

            model_data$max_std_diff <- max_std_diff
            model_data$mean_std_diff <- mean_std_diff
            model_data$n_balanced_vars <- n_balanced
            model_data$n_total_vars <- n_total

            if (verbose) {
              cat("    Balance summary (avg across imputations): Max |std diff| =", round(max_std_diff, 3),
                  ", Mean |std diff| =", round(mean_std_diff, 3), "\n")
            }
          }

        } else if ("balance" %in% names(model_result) && "Balance" %in% names(model_result$balance)) {
          # Use single imputation balance
          balance_df <- as.data.frame(model_result$balance$Balance)

          # Check for standardized difference column
          std_diff_col <- NULL
          if ("Diff.Adj" %in% names(balance_df)) {
            std_diff_col <- "Diff.Adj"
          } else if ("Diff.Target.Adj" %in% names(balance_df)) {
            std_diff_col <- "Diff.Target.Adj"
          }

          if (!is.null(std_diff_col) && nrow(balance_df) > 0) {
            # Calculate overall balance metrics
            std_diffs <- balance_df[[std_diff_col]]
            std_diffs <- std_diffs[!is.na(std_diffs)]

            if (length(std_diffs) > 0) {
              max_std_diff <- max(abs(std_diffs))
              mean_std_diff <- mean(abs(std_diffs))
              n_balanced <- sum(abs(std_diffs) < 0.1)
              n_total <- length(std_diffs)

              model_data$max_std_diff <- max_std_diff
              model_data$mean_std_diff <- mean_std_diff
              model_data$n_balanced_vars <- n_balanced
              model_data$n_total_vars <- n_total

              if (verbose) {
                cat("    Balance summary: Max |std diff| =", round(max_std_diff, 3),
                    ", Mean |std diff| =", round(mean_std_diff, 3), "\n")
                cat("    Variables balanced (<0.1):", n_balanced, "out of", n_total, "\n")
              }
            } else {
              if (verbose) cat("    No valid standardized differences found\n")
            }
          } else {
            if (verbose) cat("    Standardized difference column not found in balance table\n")
          }
        } else {
          if (verbose) cat("    No balance table available for this model\n")
        }
      }

      # Filter out intercept if requested
      if (exclude_intercept) {
        model_data <- model_data %>%
          dplyr::filter(variable != "(Intercept)")
      }

      # Filter to only the primary treatment variable
      model_data <- model_data %>%
        dplyr::filter(variable == primary_treatment)

      # Add variable type for analysis
      model_data <- model_data %>%
        dplyr::mutate(
          variable_type = "Primary Treatment"
        )

      all_results[[i]] <- model_data

    }, error = function(e) {
      if (verbose) cat("Error processing", model_name, ":", e$message, "\n")
      all_results[[i]] <- NULL
    })
  }

  # Combine all results
  combined_results <- dplyr::bind_rows(all_results)

  if (nrow(combined_results) == 0) {
    stop("No valid results found. Check your model inputs.")
  }

  # Add significance indicators for different alpha levels
  for (alpha in alpha_levels) {
    col_name <- paste0("sig_", gsub("\\.", "_", as.character(alpha)))
    combined_results[[col_name]] <- combined_results$glm_p_value < alpha
  }

  # Add bootstrap significance (CI doesn't include 0)
  if (include_bootstrap && "bootstrap_conf_low" %in% names(combined_results)) {
    combined_results$bootstrap_significant <-
      (combined_results$bootstrap_conf_low > 0 & combined_results$bootstrap_conf_high > 0) |
      (combined_results$bootstrap_conf_low < 0 & combined_results$bootstrap_conf_high < 0)
  }

  # Create summary table - focusing on primary treatment
  summary_table <- combined_results %>%
    dplyr::group_by(model_name) %>%
    dplyr::summarise(
      n_variables = dplyr::n(),
      sample_size = first(sample_size),
      across(dplyr::starts_with("sig_"), \(x) sum(x, na.rm = TRUE)),
      bootstrap_significant_count = if("bootstrap_significant" %in% names(.)) sum(bootstrap_significant, na.rm = TRUE) else NA,
      min_p_value = min(glm_p_value, na.rm = TRUE),
      # Focus on primary treatment significance
      primary_treatment_significant = any(glm_p_value < 0.05, na.rm = TRUE),
      primary_treatment_p_value = min(glm_p_value, na.rm = TRUE),
      primary_treatment_estimate = first(glm_estimate),
      # Add model fit statistics
      AIC = first(AIC),
      BIC = first(BIC),
      BICc = first(BICc),
      # Add MI info
      n_imputations = first(n_imputations),
      inference_method = first(inference_method),
      .groups = "drop"
    )

  # Filter results based on focus_on_treatment_only setting
  if (focus_on_treatment_only) {
    # Only include models where primary treatment is significant
    sig_details <- combined_results %>%
      dplyr::filter(glm_p_value < 0.05) %>%
      dplyr::arrange(glm_p_value)

    if (verbose) {
      significant_models <- length(unique(sig_details$model_name))
      total_models <- length(unique(combined_results$model_name))
      cat("Found", significant_models, "model(s) with significant primary treatment effects (p < 0.05) out of", total_models, "total models\n")
      if (significant_models > 0) {
        cat("Models with significant treatment:", paste(unique(sig_details$model_name), collapse = ", "), "\n")
      }
    }
  } else {
    # Include all treatment results regardless of significance
    sig_details <- combined_results %>%
      dplyr::arrange(glm_p_value)

    if (verbose) {
      cat("Showing primary treatment results from all", nrow(combined_results), "models\n")
    }
  }

  # Select columns for output
  base_columns <- c("model_name", "variable", "glm_estimate", "glm_se", "glm_p_value")

  # Add confidence intervals if they exist
  if (all(c("glm_conf_low", "glm_conf_high") %in% names(sig_details))) {
    base_columns <- c(base_columns, "glm_conf_low", "glm_conf_high")
  }

  # Add bootstrap columns if they exist
  bootstrap_columns <- c("bootstrap_estimate", "bootstrap_conf_low", "bootstrap_conf_high")
  existing_bootstrap <- bootstrap_columns[bootstrap_columns %in% names(sig_details)]
  base_columns <- c(base_columns, existing_bootstrap)

  # Add sample size if it exists
  if ("sample_size" %in% names(sig_details)) {
    base_columns <- c(base_columns, "sample_size")
  }

  # Add MI-specific columns if they exist
  mi_columns <- c("n_imputations", "inference_method", "fmi", "df_adjusted")
  existing_mi <- mi_columns[mi_columns %in% names(sig_details)]
  base_columns <- c(base_columns, existing_mi)

  # Add balance statistics if they exist
  balance_columns <- c("max_std_diff", "mean_std_diff", "n_balanced_vars", "n_total_vars")
  existing_balance <- balance_columns[balance_columns %in% names(sig_details)]
  base_columns <- c(base_columns, existing_balance)

  # Add model fit statistics if they exist
  fit_columns <- c("AIC", "BIC", "BICc")
  existing_fit <- fit_columns[fit_columns %in% names(sig_details)]
  base_columns <- c(base_columns, existing_fit)

  # Add bootstrap significance if it exists and requested
  if(include_bootstrap && "bootstrap_significant" %in% names(sig_details)) {
    base_columns <- c(base_columns, "bootstrap_significant")
  }

  # Add significance flags
  sig_columns <- names(sig_details)[grepl("^sig_", names(sig_details))]
  base_columns <- c(base_columns, sig_columns)

  # Select only existing columns
  sig_details <- sig_details %>%
    dplyr::select(dplyr::all_of(base_columns))

  # Print summary if verbose
  if (verbose) {
    cat("\n=== TREATMENT SIGNIFICANCE SUMMARY ===\n")
    cat(paste("Primary treatment variable:", primary_treatment, "\n"))

    # Show inference method being used
    inference_methods <- unique(combined_results$inference_method)
    cat(paste("Inference method:", paste(inference_methods, collapse = ", "), "\n\n"))

    # Treatment significance summary
    cat("Treatment Significance Across Models:\n")

    treatment_summary <- summary_table %>%
      dplyr::arrange(primary_treatment_p_value)

    for (i in 1:nrow(treatment_summary)) {
      status <- if (treatment_summary$primary_treatment_significant[i]) "yes" else "no"
      mi_info <- if (treatment_summary$n_imputations[i] > 1)
        paste0(" (", treatment_summary$n_imputations[i], " imputations)") else ""

      cat(status, " ", treatment_summary$model_name[i], mi_info, ": p = ",
          round(treatment_summary$primary_treatment_p_value[i], 4),
          ", effect = ", round(treatment_summary$primary_treatment_estimate[i], 4), "\n")
    }

    significant_count <- sum(treatment_summary$primary_treatment_significant)
    cat("\nTotal models with significant treatment effect:", significant_count,
        "out of", nrow(treatment_summary), "\n")

    # Model fit comparison if available
    if ("AIC" %in% names(summary_table) && any(!is.na(summary_table$AIC))) {
      cat("\nModel Fit Statistics:\n")
      fit_summary <- summary_table %>%
        dplyr::select(model_name, AIC, BIC, BICc) %>%
        dplyr::arrange(AIC)

      for (i in 1:nrow(fit_summary)) {
        cat("  ", fit_summary$model_name[i], ": AIC =", round(fit_summary$AIC[i], 2),
            ", BIC =", round(fit_summary$BIC[i], 2),
            ", BICc =", round(fit_summary$BICc[i], 2), "\n")
      }

      cat("\n  Best fitting model (lowest AIC):", fit_summary$model_name[1], "\n")
    }

    # Balance assessment summary if available
    if ("max_std_diff" %in% names(combined_results) && any(!is.na(combined_results$max_std_diff))) {
      cat("\nBalance Assessment Summary:\n")

      balance_summary <- combined_results %>%
        dplyr::select(model_name, max_std_diff, mean_std_diff, n_balanced_vars, n_total_vars) %>%
        dplyr::distinct() %>%
        dplyr::filter(!is.na(max_std_diff)) %>%
        dplyr::arrange(max_std_diff)

      if (nrow(balance_summary) > 0) {
        for (i in 1:nrow(balance_summary)) {
          balance_status <- if (balance_summary$max_std_diff[i] < 0.1) {
            " Good"
          } else if (balance_summary$max_std_diff[i] < 0.25) {
            " Moderate"
          } else {
            "Warning: Poor"
          }

          cat("  ", balance_summary$model_name[i], ": Max |std diff| =",
              round(balance_summary$max_std_diff[i], 3),
              ", Balanced vars:", balance_summary$n_balanced_vars[i], "/", balance_summary$n_total_vars[i],
              balance_status, "\n")
        }
      }
    }

    # Best result summary
    if (significant_count > 0) {
      best_result <- treatment_summary %>%
        dplyr::filter(primary_treatment_significant == TRUE) %>%
        dplyr::slice_min(primary_treatment_p_value, n = 1)

      cat("\nMost significant treatment result:\n")
      cat("  Model:", best_result$model_name[1], "\n")
      cat("  P-value:", round(best_result$primary_treatment_p_value[1], 4), "\n")
      cat("  Effect size:", round(best_result$primary_treatment_estimate[1], 4), "\n")

      # Show if using MI
      if (best_result$n_imputations[1] > 1) {
        cat("  Multiple imputations:", best_result$n_imputations[1], "\n")
      }
    }
  }

  # Return comprehensive results
  return(list(
    detailed_results = combined_results,
    summary_table = summary_table,
    significant_only = sig_details,
    alpha_levels = alpha_levels,
    n_models = length(model_results_list),
    model_names = model_names
  ))
}


#' Quick Significance Check
quick_sig_check <- function(model_results_list, model_names = NULL, alpha = 0.05,
                            primary_treatment = "treatment", focus_on_treatment_only = TRUE,
                            use_pooled = TRUE) {

  results <- check_treatment_significance(
    model_results_list = model_results_list,
    model_names = model_names,
    primary_treatment = primary_treatment,
    alpha_levels = alpha,
    focus_on_treatment_only = focus_on_treatment_only,
    use_pooled = use_pooled,
    verbose = FALSE
  )

  if (focus_on_treatment_only) {
    sig_models <- results$summary_table %>%
      dplyr::filter(primary_treatment_significant == TRUE)

    if (nrow(sig_models) > 0) {
      cat("Models with significant", primary_treatment, "effects (p <", alpha, "):\n")
      print(sig_models %>% dplyr::select(model_name, primary_treatment_p_value, primary_treatment_estimate, sample_size, n_imputations, AIC, BIC, BICc))
    } else {
      cat("No models with significant", primary_treatment, "effects at p <", alpha, "\n")
    }
  } else {
    # Show all treatment results
    all_models <- results$summary_table
    cat("Primary treatment results for all models:\n")
    print(all_models %>% dplyr::select(model_name, primary_treatment_p_value, primary_treatment_estimate, sample_size, n_imputations, AIC, BIC, BICc))
  }

  return(results$summary_table)
}
