#' Quick Model Coefficient Table
#'
#' A simple helper function to quickly view coefficients (b-values) and p-values
#' from a fitted model during interactive analysis. Designed for quick inspection,
#' not publication-ready tables.
#'
#' @param model A fitted model object (glm, lm, or similar) OR a list from cbps_weighted_analysis/npcbps_weighted_analysis
#' @param digits Integer. Number of decimal places to round to (default: 3)
#' @param include_intercept Logical. Whether to include the intercept term (default: FALSE)
#' @param sort_by Character. Sort by "none", "pvalue", "b_value", or "abs_b_value" (default: "none")
#' @param sig_only Logical. Show only significant results (default: FALSE)
#' @param alpha Numeric. Significance level if sig_only = TRUE (default: 0.05)
#'
#' @return A data frame with columns: variable, b_value, p_value (and optionally sig)
#'
#' @details
#' This is a lightweight function for quick model inspection during interactive
#' analysis. For comprehensive results extraction and publication tables, use
#' \code{\link{extract_cbps_results}} instead.
#'
#' @examples
#' \dontrun{
#' # Quick look at any model
#' model <- glm(outcome ~ treatment + age, data = mydata)
#' quick_coef_table(model)
#'
#' # From CBPS analysis results
#' results <- cbps_weighted_analysis(...)
#' quick_coef_table(results)
#' quick_coef_table(results$model)  # same thing
#'
#' # Only significant results, sorted by p-value
#' quick_coef_table(model, sig_only = TRUE, sort_by = "pvalue")
#'
#' # Sort by effect size
#' quick_coef_table(model, sort_by = "abs_b_value", include_intercept = FALSE)
#' }
#'
#' @seealso \code{\link{extract_cbps_results}} for comprehensive results extraction
#' @export
quick_coef_table <- function(model,
                             digits = 3,
                             include_intercept = FALSE,
                             sort_by = "none",
                             sig_only = FALSE,
                             alpha = 0.05) {

  # Check if input is from extract_cbps_results()
  if (is.list(model) && "estimates" %in% names(model)) {
    # Input from extract_cbps_results() - use estimates dataframe
    estimates_df <- model$estimates

    # Create result dataframe from estimates
    result <- data.frame(
      variable = estimates_df$variable,
      b_value = round(estimates_df$glm_estimate, digits),
      p_value = round(estimates_df$glm_p_value, digits),
      stringsAsFactors = FALSE
    )

  } else if (is.list(model) && "model" %in% names(model)) {
    # Input from cbps/npcbps analysis functions - extract model
    model <- model$model

    # Extract coefficient summary
    coef_summary <- summary(model)$coefficients

    # Create simple data frame
    result <- data.frame(
      variable = rownames(coef_summary),
      b_value = round(coef_summary[, "Estimate"], digits),
      p_value = round(coef_summary[, "Pr(>|t|)"], digits),
      stringsAsFactors = FALSE
    )

  } else {
    # Direct model object
    # Extract coefficient summary
    coef_summary <- summary(model)$coefficients

    # Create simple data frame
    result <- data.frame(
      variable = rownames(coef_summary),
      b_value = round(coef_summary[, "Estimate"], digits),
      p_value = round(coef_summary[, "Pr(>|t|)"], digits),
      stringsAsFactors = FALSE
    )
  }

  # Add significance markers
  result$sig <- ifelse(result$p_value < 0.001, "***",
                       ifelse(result$p_value < 0.01, "**",
                              ifelse(result$p_value < 0.05, "*", "")))

  # Filter for significant only if requested
  if (sig_only) {
    result <- result[result$p_value < alpha, ]
  }

  # Remove intercept if requested
  if (!include_intercept) {
    result <- result[result$variable != "(Intercept)", ]
  }

  # Sort if requested
  if (sort_by == "pvalue") {
    result <- result[order(result$p_value), ]
  } else if (sort_by == "b_value") {
    result <- result[order(result$b_value), ]
  } else if (sort_by == "abs_b_value") {
    result <- result[order(-abs(result$b_value)), ]
  }

  # Reset row names for cleaner display
  rownames(result) <- NULL

  return(result)
}


#' Print Quick Coefficient Table
#'
#' Prints the quick coefficient table with nice formatting to console.
#' This is a convenience wrapper around \code{quick_coef_table} that prints
#' the results with a header for easier reading during interactive analysis.
#'
#' @param model A fitted model object (glm, lm, or similar) OR a list from cbps_weighted_analysis OR output from extract_cbps_results
#' @param digits Integer. Number of decimal places to round to (default: 3)
#' @param include_intercept Logical. Whether to include the intercept term (default: FALSE)
#' @param sort_by Character. Sort by "none", "pvalue", "b_value", or "abs_b_value" (default: "none")
#' @param sig_only Logical. Show only significant results (default: FALSE)
#' @param alpha Numeric. Significance level if sig_only = TRUE (default: 0.05)
#'
#' @return Invisibly returns the data frame (mainly used for side effect of printing)
#'
#' @examples
#' \dontrun{
#' results <- cbps_weighted_analysis(...)
#' print_coef_table(results)
#' print_coef_table(results, sig_only = TRUE)
#'
#' # Or with extracted results
#' extracted <- extract_cbps_results(results)
#' print_coef_table(extracted)
#' }
#'
#' @export
print_coef_table <- function(model,
                             digits = 3,
                             include_intercept = FALSE,
                             sort_by = "none",
                             sig_only = FALSE,
                             alpha = 0.05) {

  result <- quick_coef_table(model, digits, include_intercept, sort_by, sig_only, alpha)

  cat("\n")
  cat("═══════════════════════════════════════\n")
  cat("  Quick Model Coefficient Summary\n")
  cat("═══════════════════════════════════════\n")
  if (sig_only) {
    cat(sprintf("  Showing only p < %.2f\n", alpha))
  }
  cat("  Significance: *** p<0.001, ** p<0.01, * p<0.05\n")
  cat("───────────────────────────────────────\n\n")

  print(result, row.names = FALSE)

  cat("\n───────────────────────────────────────\n")
  cat(sprintf("  Total variables: %d\n", nrow(result)))
  if (!sig_only) {
    n_sig <- sum(result$sig != "")
    cat(sprintf("  Significant (p<0.05): %d\n", n_sig))
  }
  cat("═══════════════════════════════════════\n\n")

  invisible(result)
}


#' Ultra-Quick Model Check (Alias)
#'
#' Ultra-short alias for quick_coef_table for even faster interactive use.
#' Just type \code{qm(model)} to quickly check your model coefficients.
#'
#' @param model A fitted model object or CBPS results list or extract_cbps_results output
#' @param ... Additional arguments passed to \code{quick_coef_table}
#'
#' @return A data frame with variable names, coefficients, and p-values
#'
#' @examples
#' \dontrun{
#' # Super quick inspection
#' qm(model)
#' qm(cbps_results)
#' qm(extract_cbps_results(cbps_results))
#'
#' # With options
#' qm(model, sig_only = TRUE)
#' qm(model, sort_by = "pvalue")
#' }
#'
#' @export
qm <- function(model, ...) {
  quick_coef_table(model, ...)
}
