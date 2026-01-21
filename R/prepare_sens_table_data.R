#' Prepare Sensitivity Analysis Table Data
#'
#' Converts cbps_weighted_analysis() outputs into format for create_sens_table()
#'
#' @param model_list Named list of lists. Structure:
#'   list(
#'     main = list(outcome1 = model1, outcome2 = model2, ...),
#'     sens1 = list(outcome1 = model1, outcome2 = model2, ...),
#'     ...
#'   )
#' @param sens_labels Named character vector. Labels for sensitivity types.
#'   Example: c(main = "Main Model", sens1 = "Adjusted for Z")
#' @param treatment_var Character. Name of treatment variable to extract from models.
#'   If NULL, uses first non-intercept coefficient.
#' @param verbose Logical. Print progress messages (default: TRUE)
#'
#' @return Data frame formatted for create_sens_table()
#'
#' @examples
#' \dontrun{
#' # Run models
#' main_outcome1 <- cbps_weighted_analysis(...)
#' main_outcome2 <- cbps_weighted_analysis(...)
#' sens_outcome1 <- cbps_weighted_analysis(..., additional_predictors = c("z"))
#' sens_outcome2 <- cbps_weighted_analysis(..., additional_predictors = c("z"))
#'
#' # Prepare data
#' sens_data <- prepare_sens_table_data(
#'   model_list = list(
#'     main = list(outcome1 = main_outcome1, outcome2 = main_outcome2),
#'     adj_z = list(outcome1 = sens_outcome1, outcome2 = sens_outcome2)
#'   ),
#'   sens_labels = c(main = "Main Model", adj_z = "Adjusted for Covariate Z"),
#'   treatment_var = "pregnancy_advice"
#' )
#'
#' # Create table
#' create_sens_table(sens_data)
#' }
#'
#' @export
prepare_sens_table_data <- function(
    model_list,
    sens_labels = NULL,
    treatment_var = NULL,
    use_pooled = TRUE,
    verbose = TRUE
) {

  required_packages <- c("dplyr", "purrr", "tibble")
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  }

  # Check that model_list is structured correctly
  if (!is.list(model_list) || is.null(names(model_list))) {
    stop("model_list must be a named list of lists.\n",
         "Example structure:\n",
         "list(\n",
         "  main = list(outcome1 = model1, outcome2 = model2),\n",
         "  sens1 = list(outcome1 = model1, outcome2 = model2)\n",
         ")")
  }

  # Set default sensitivity labels if not provided
  if (is.null(sens_labels)) {
    sens_labels <- setNames(names(model_list), names(model_list))
  }

  # Helper function to extract info from a single model
  extract_model_info <- function(model_result, sens_type, outcome_name, treat_var) {

    tryCatch({

      # Check if pooled results exist and should be used
      use_pooled_for_model <- use_pooled && !is.null(model_result$pooled_results)

      if (use_pooled_for_model) {
        # USE POOLED RESULTS (Rubin's rules)
        if (verbose) cat("      Using Rubin's rules pooled results\n")

        pooled_data <- model_result$pooled_results

        # If treatment_var not specified, use the first non-intercept term
        if (is.null(treat_var)) {
          available_terms <- pooled_data$term[pooled_data$term != "(Intercept)"]
          if (length(available_terms) == 0) {
            stop("No non-intercept coefficients found in pooled results")
          }
          treat_var <- available_terms[1]
          if (verbose) {
            cat("      Using", treat_var, "as treatment variable (first non-intercept term)\n")
          }
        }

        # Check if treatment variable exists
        if (!treat_var %in% pooled_data$term) {
          stop(paste("Treatment variable", treat_var, "not found in pooled results for", outcome_name,
                     "\nAvailable coefficients:", paste(pooled_data$term, collapse = ", ")))
        }

        # Extract pooled info
        pooled_row <- pooled_data %>% filter(term == treat_var)

        glm_estimate <- pooled_row$estimate
        glm_se <- pooled_row$se
        glm_p_value <- pooled_row$p.value
        glm_conf_low <- pooled_row$ci_lower
        glm_conf_high <- pooled_row$ci_upper
        fmi <- pooled_row$fmi
        df_adjusted <- pooled_row$df
        inference_method <- "Rubin's Rules (MI)"
        n_imputations <- model_result$n_imputations

      } else {
        # USE SINGLE IMPUTATION (backward compatible)
        if (verbose && !is.null(model_result$pooled_results)) {
          cat("      Using single imputation results (use_pooled=FALSE)\n")
        }

        # Get GLM summary
        glm_summary <- summary(model_result$model)
        coef_table <- glm_summary$coefficients

        # If treatment_var not specified, use the first non-intercept term
        if (is.null(treat_var)) {
          available_terms <- rownames(coef_table)[rownames(coef_table) != "(Intercept)"]
          if (length(available_terms) == 0) {
            stop("No non-intercept coefficients found in model")
          }
          treat_var <- available_terms[1]
          if (verbose) {
            cat("      Using", treat_var, "as treatment variable (first non-intercept term)\n")
          }
        }

        # Check if treatment variable exists in model
        if (!treat_var %in% rownames(coef_table)) {
          stop(paste("Treatment variable", treat_var, "not found in model for", outcome_name,
                     "\nAvailable coefficients:", paste(rownames(coef_table), collapse = ", ")))
        }

        # Extract GLM info
        glm_estimate <- coef_table[treat_var, "Estimate"]
        glm_se <- coef_table[treat_var, "Std. Error"]
        glm_p_value <- coef_table[treat_var, "Pr(>|t|)"]

        # Calculate GLM confidence intervals
        glm_conf_low <- glm_estimate - 1.96 * glm_se
        glm_conf_high <- glm_estimate + 1.96 * glm_se

        fmi <- NA_real_
        df_adjusted <- NA_real_
        inference_method <- "Single Imputation"
        n_imputations <- 1
      }

      # Extract bootstrap info if available
      bootstrap_available <- !is.null(model_result$bootstrap_summary) &&
        nrow(model_result$bootstrap_summary) > 0

      if (bootstrap_available) {
        boot_row <- model_result$bootstrap_summary %>%
          filter(term == treat_var)

        if (nrow(boot_row) > 0) {
          bootstrap_estimate <- boot_row$estimate
          bootstrap_conf_low <- boot_row$ci_lower
          bootstrap_conf_high <- boot_row$ci_upper
        } else {
          bootstrap_estimate <- NA_real_
          bootstrap_conf_low <- NA_real_
          bootstrap_conf_high <- NA_real_
        }
      } else {
        bootstrap_estimate <- NA_real_
        bootstrap_conf_low <- NA_real_
        bootstrap_conf_high <- NA_real_
      }

      # Get sample size
      sample_size <- model_result$sample_sizes$final

      # Return as tibble
      tibble(
        sensitivity_type = sens_type,
        model_name = outcome_name,
        variable = treat_var,
        glm_estimate = glm_estimate,
        glm_se = glm_se,
        glm_p_value = glm_p_value,
        glm_conf_low = glm_conf_low,
        glm_conf_high = glm_conf_high,
        bootstrap_estimate = bootstrap_estimate,
        bootstrap_conf_low = bootstrap_conf_low,
        bootstrap_conf_high = bootstrap_conf_high,
        sample_size = sample_size,
        n_imputations = n_imputations,
        inference_method = inference_method,
        fmi = fmi,
        df_adjusted = df_adjusted
      )

    }, error = function(e) {
      warning("Error extracting info for ", outcome_name, " in ", sens_type, ": ", e$message)
      return(NULL)
    })
  }

  # Process all models
  if (verbose) cat("Extracting information from models...\n")

  all_data <- map_dfr(names(model_list), function(sens_type) {

    if (verbose) cat("  Processing sensitivity type:", sens_type, "\n")

    sens_models <- model_list[[sens_type]]

    # Check that this is a list of models
    if (!is.list(sens_models)) {
      warning("Sensitivity type ", sens_type, " is not a list. Skipping.")
      return(NULL)
    }

    # Check for names
    if (is.null(names(sens_models))) {
      stop("Models within sensitivity type ", sens_type, " must be named.\n",
           "Example: list(outcome1 = model1, outcome2 = model2)")
    }

    # Extract info from each outcome model
    map_dfr(names(sens_models), function(outcome_name) {
      if (verbose) cat("    Extracting:", outcome_name, "\n")
      extract_model_info(
        model_result = sens_models[[outcome_name]],
        sens_type = sens_type,
        outcome_name = outcome_name,
        treat_var = treatment_var
      )
    })
  })

  # Check if any data was extracted
  if (is.null(all_data) || nrow(all_data) == 0) {
    stop("No data could be extracted from models. Check model structure and error messages.")
  }

  if (verbose) {
    cat("\nSuccessfully extracted data for:\n")
    cat("  Sensitivity types:", length(unique(all_data$sensitivity_type)), "\n")
    cat("  Outcomes:", length(unique(all_data$model_name)), "\n")
    cat("  Total rows:", nrow(all_data), "\n")

    # Show inference method info
    inference_methods <- unique(all_data$inference_method)
    cat("  Inference methods:", paste(inference_methods, collapse = ", "), "\n")

    mi_models <- sum(all_data$n_imputations > 1)
    if (mi_models > 0) {
      cat("  Models using multiple imputation:", mi_models, "\n")
    }
  }

  return(all_data)
}
