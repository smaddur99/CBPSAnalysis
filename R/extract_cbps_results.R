#' Extract CBPS Results into Structured Format
#'
#' Extracts key results from CBPS analysis into a structured format suitable for
#' visualization and further analysis.
#'
#' @param cbps_results List. Output from \code{\link{cbps_weighted_analysis}}
#' @param include_balance_details Logical. Whether to include detailed balance statistics (default: TRUE)
#'
#' @return If \code{include_balance_details = TRUE}, returns a list containing:
#' \itemize{
#'   \item \code{estimates}: Combined GLM and bootstrap estimates
#'   \item \code{balance}: Detailed balance statistics
#'   \item \code{sample_sizes}: Sample size information
#'   \item \code{summary_stats}: Summary balance metrics
#' }
#'
#' If \code{include_balance_details = FALSE}, returns a data frame with model estimates and basic information.
#'
#' @details
#' This function processes the complex output from CBPS analysis and creates clean,
#' structured data frames suitable for:
#' \itemize{
#'   \item Creating visualization plots
#'   \item Comparing multiple models
#'   \item Generating publication tables
#'   \item Further statistical analysis
#' }
#'
#' The function automatically handles different column naming conventions from
#' different versions of the WeightIt and cobalt packages.
#'
#' @examples
#' \dontrun{
#' # After running CBPS analysis
#' results <- cbps_weighted_analysis(...)
#'
#' # Extract detailed results
#' detailed <- extract_cbps_results(results, include_balance_details = TRUE)
#'
#' # Access components
#' detailed$estimates        # Model and bootstrap estimates
#' detailed$balance         # Balance statistics
#' detailed$sample_sizes    # Sample size info
#'
#' # Extract simple results for visualization
#' simple <- extract_cbps_results(results, include_balance_details = FALSE)
#' }
#'
#' @seealso \code{\link{cbps_weighted_analysis}}, \code{\link{quick_results_summary}}
#' @export
#' @importFrom dplyr select all_of filter mutate full_join rename
#' @importFrom tibble tibble as_tibble rownames_to_column
#' @importFrom broom tidy
extract_cbps_results <- function(cbps_results, include_balance_details = TRUE) {

  # Load required libraries
  required_packages <- c("dplyr", "tibble", "broom")
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  }

  # 1. EXTRACT MODEL ESTIMATES
  model_summary <- broom::tidy(cbps_results$model, conf.int = TRUE) %>%
    mutate(
      analysis_type = "GLM_Model",
      sample_size = cbps_results$sample_sizes$final
    ) %>%
    rename(
      variable = term,
      glm_estimate = estimate,
      glm_se = std.error,
      glm_statistic = statistic,
      glm_p_value = p.value,
      glm_conf_low = conf.low,
      glm_conf_high = conf.high
    )

  # 2. EXTRACT BOOTSTRAP RESULTS
  bootstrap_results <- cbps_results$bootstrap_summary %>%
    mutate(
      analysis_type = "Bootstrap",
      sample_size = cbps_results$sample_sizes$final
    ) %>%
    rename(
      variable = term,
      bootstrap_estimate = estimate,
      bootstrap_se = se,
      bootstrap_conf_low = ci_lower,
      bootstrap_conf_high = ci_upper
    )

  # 3. EXTRACT BALANCE STATISTICS
  balance_data <- cbps_results$balance

  # Extract balance table as dataframe
  if ("Balance" %in% names(balance_data)) {
    balance_df <- balance_data$Balance %>%
      as.data.frame() %>%
      rownames_to_column("variable") %>%
      as_tibble() %>%
      mutate(
        analysis_type = "Balance",
        sample_size = cbps_results$sample_sizes$final
      )

    # Clean up column names (they vary depending on WeightIt version)
    balance_names <- names(balance_df)
    new_names <- balance_names

    # Standard renaming patterns for balance table structure
    if ("Diff.Un" %in% balance_names) new_names[new_names == "Diff.Un"] <- "diff_unadj"
    if ("Diff.Adj" %in% balance_names) new_names[new_names == "Diff.Adj"] <- "diff_adj"
    if ("Diff.Target.Adj" %in% balance_names) new_names[new_names == "Diff.Target.Adj"] <- "diff_adj"
    if ("Corr.Adj" %in% balance_names) new_names[new_names == "Corr.Adj"] <- "corr_adj"
    if ("M.Threshold" %in% balance_names) new_names[new_names == "M.Threshold"] <- "balance_status"
    if ("M.0.Un" %in% balance_names) new_names[new_names == "M.0.Un"] <- "mean_control_unadj"
    if ("M.1.Un" %in% balance_names) new_names[new_names == "M.1.Un"] <- "mean_treated_unadj"
    if ("M.0.Adj" %in% balance_names) new_names[new_names == "M.0.Adj"] <- "mean_control_adj"
    if ("M.1.Adj" %in% balance_names) new_names[new_names == "M.1.Adj"] <- "mean_treated_adj"
    if ("V.Ratio.Un" %in% balance_names) new_names[new_names == "V.Ratio.Un"] <- "var_ratio_unadj"
    if ("V.Ratio.Adj" %in% balance_names) new_names[new_names == "V.Ratio.Adj"] <- "var_ratio_adj"

    names(balance_df) <- new_names

  } else {
    # Fallback if balance structure is different
    balance_df <- tibble(
      variable = "Balance data not available",
      analysis_type = "Balance",
      sample_size = cbps_results$sample_sizes$final
    )
  }

  # 4. COMBINE ALL RESULTS
  # Join model and bootstrap results
  combined_estimates <- model_summary %>%
    full_join(bootstrap_results, by = c("variable", "sample_size")) %>%
    mutate(analysis_type = "Combined")

  # 5. CREATE FINAL RESULTS DATAFRAME
  if (include_balance_details) {
    # Include detailed balance statistics
    final_results <- list(
      estimates = combined_estimates,
      balance = balance_df,
      sample_sizes = tibble(
        metric = c("initial_n", "final_n", "n_dropped"),
        value = c(
          cbps_results$sample_sizes$initial,
          cbps_results$sample_sizes$final,
          cbps_results$sample_sizes$initial - cbps_results$sample_sizes$final
        )
      ),
      # Summary statistics
      summary_stats = tibble(
        metric = c(
          "max_abs_std_diff_before",
          "max_abs_std_diff_after",
          "n_variables_balanced",
          "mean_absolute_std_diff_before",
          "mean_absolute_std_diff_after"
        ),
        value = c(
          if ("diff_unadj" %in% names(balance_df)) max(abs(balance_df$diff_unadj), na.rm = TRUE) else NA,
          if ("diff_adj" %in% names(balance_df)) max(abs(balance_df$diff_adj), na.rm = TRUE) else NA,
          if ("diff_adj" %in% names(balance_df)) sum(abs(balance_df$diff_adj) < 0.1, na.rm = TRUE) else NA,
          if ("diff_unadj" %in% names(balance_df)) mean(abs(balance_df$diff_unadj), na.rm = TRUE) else NA,
          if ("diff_adj" %in% names(balance_df)) mean(abs(balance_df$diff_adj), na.rm = TRUE) else NA
        )
      )
    )
  } else {
    # Simple version - just estimates and basic info
    final_results <- combined_estimates %>%
      dplyr::mutate(
        initial_sample_size = cbps_results$sample_sizes$initial,
        final_sample_size = cbps_results$sample_sizes$final,
        n_dropped = cbps_results$sample_sizes$initial - cbps_results$sample_sizes$final
      )
  }

  return(final_results)
}


#' Quick Results Summary for Visualization
#'
#' Creates a simplified summary of CBPS results optimized for visualization.
#'
#' @param cbps_results List. Output from \code{\link{cbps_weighted_analysis}}
#'
#' @return A data frame with key estimates and confidence intervals, ready for plotting
#'
#' @details
#' This is a convenience function that extracts only the essential columns needed
#' for creating plots and quick summaries. It's equivalent to calling
#' \code{extract_cbps_results()} with \code{include_balance_details = FALSE} and
#' selecting key columns.
#'
#' @examples
#' \dontrun{
#' # Quick visualization-ready data
#' viz_data <- quick_results_summary(results)
#'
#' # Create a simple plot
#' library(ggplot2)
#' viz_data %>%
#'   filter(variable != "(Intercept)") %>%
#'   ggplot(aes(x = variable, y = glm_estimate)) +
#'   geom_point() +
#'   geom_errorbar(aes(ymin = glm_conf_low, ymax = glm_conf_high))
#' }
#'
#' @seealso \code{\link{extract_cbps_results}}
#' @export
quick_results_summary <- function(cbps_results) {
  result <- extract_cbps_results(cbps_results, include_balance_details = FALSE)

  # Debug: show what columns we actually have
  cat("Available columns:", paste(names(result), collapse = ", "), "\n")

  # Select only columns that exist
  available_cols <- names(result)
  desired_cols <- c(
    "variable",
    "glm_estimate", "glm_se", "glm_p_value", "glm_conf_low", "glm_conf_high",
    "bootstrap_estimate", "bootstrap_se", "bootstrap_conf_low", "bootstrap_conf_high",
    "final_sample_size"
  )

  # Only select columns that exist
  cols_to_select <- desired_cols[desired_cols %in% available_cols]

  return(result %>% dplyr::select(all_of(cols_to_select)))
}
