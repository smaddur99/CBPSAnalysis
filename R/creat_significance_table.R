#' Create Publication-Ready Significance Table
#'
#' Creates a formatted table of treatment effects across model specifications,
#' suitable for publication or presentation.
#'
#' @param sig_results List or data frame. Output from \code{\link{check_treatment_significance}} or equivalent data frame
#' @param table_title Character. Title for the table (default: "Treatment Effects Across Model Specifications")
#' @param filename Character. Optional filename to save the table (default: NULL)
#' @param treatment_label Character. Label for the treatment variable (default: NULL, uses variable name)
#' @param model_labels Named character vector. Pretty names for models (default: NULL)
#' @param include_bootstrap Logical. Whether to include bootstrap results (default: TRUE)
#' @param decimal_places Integer. Number of decimal places for estimates (default: 3)
#' @param font_size Integer. Font size for the table (default: 16)
#' @param font_family Character. Font family for the table (default: "Times New Roman")
#' @param show_all_models Logical. Whether to show all models or just significant ones (default: FALSE)
#' @param alpha_threshold Numeric. Threshold for highlighting significance (default: 0.05)
#'
#' @return A gt table object that can be displayed or saved
#'
#' @details
#' This function creates publication-ready tables using the gt package with:
#' \itemize{
#'   \item Professional formatting suitable for academic papers
#'   \item Significance highlighting (green background for significant results)
#'   \item Confidence intervals for both GLM and bootstrap estimates
#'   \item Balance assessment ratios
#'   \item Sample size information
#'   \item Customizable styling options
#' }
#'
#' The function automatically handles:
#' \itemize{
#'   \item Missing confidence intervals
#'   \item Bootstrap vs. non-bootstrap results
#'   \item Different balance table formats
#'   \item Significance stars (* p < 0.05, ** p < 0.01, *** p < 0.001)
#' }
#'
#' @examples
#' \dontrun{
#' # After running significance analysis
#' sig_results <- check_treatment_significance(...)
#'
#' # Create basic table
#' table1 <- create_significance_table(
#'   sig_results = sig_results,
#'   treatment_label = "Treatment Effect"
#' )
#'
#' # Create customized table with all models
#' table2 <- create_significance_table(
#'   sig_results = sig_results,
#'   table_title = "Comparison of Model Specifications",
#'   model_labels = c("Model_1" = "Basic Model", "Model_2" = "Full Model"),
#'   show_all_models = TRUE,
#'   filename = "treatment_table.html"
#' )
#' }
#'
#' @seealso \code{\link{check_treatment_significance}}, \code{\link{quick_treatment_table}}, \code{\link{full_treatment_table}}
#' @export
#' @importFrom gt gt cols_label fmt_number tab_header tab_style cell_text cell_fill cells_body cells_column_labels cells_title cols_hide cols_align tab_footnote tab_options gtsave md px
#' @importFrom dplyr select all_of mutate case_when arrange filter
create_significance_table <- function(
    sig_results,
    table_title = "Treatment Effects Across Model Specifications",
    filename = NULL,
    treatment_label = NULL,
    model_labels = NULL,
    include_bootstrap = TRUE,
    decimal_places = 3,
    font_size = 16,
    font_family = "Times New Roman",
    show_all_models = FALSE,
    alpha_threshold = 0.05
) {

  # Load required libraries
  required_packages <- c("dplyr", "gt", "tibble")
  for (pkg in required_packages) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      stop(paste("Package", pkg, "is required but not installed."))
    }
  }

  # Extract results - handle both direct dataframe and list input
  if ("significant_only" %in% names(sig_results)) {
    if (show_all_models) {
      table_data <- sig_results$detailed_results
    } else {
      table_data <- sig_results$significant_only
    }
  } else {
    table_data <- sig_results
  }

  # Debug: Check what columns are available
  cat("Available columns in table_data:", paste(names(table_data), collapse = ", "), "\n")

  # Check for required columns
  required_cols <- c("variable", "glm_estimate", "glm_p_value", "model_name")
  missing_cols <- required_cols[!required_cols %in% names(table_data)]

  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # Create pretty treatment label
  if (is.null(treatment_label)) {
    treatment_label <- unique(table_data$variable)[1]  # Use the first (should be only) variable name
  }

  # Create pretty model names
  if (is.null(model_labels)) {
    model_labels <- setNames(unique(table_data$model_name), unique(table_data$model_name))
  }

  # Prepare the data
  formatted_data <- table_data %>%
    dplyr::mutate(
      # Create pretty model names
      pretty_model = ifelse(model_name %in% names(model_labels),
                            model_labels[model_name],
                            model_name),

      # Format confidence intervals (with safety check)
      glm_conf_int = if(all(c("glm_conf_low", "glm_conf_high") %in% names(.))) {
        paste0("(",
               sprintf(paste0("%.", decimal_places, "f"), glm_conf_low),
               ", ",
               sprintf(paste0("%.", decimal_places, "f"), glm_conf_high),
               ")")
      } else {
        "CI not available"
      },

      # Format bootstrap confidence intervals if available
      bootstrap_conf_int = if(include_bootstrap && all(c("bootstrap_conf_low", "bootstrap_conf_high") %in% names(.))) {
        paste0("(",
               sprintf(paste0("%.", decimal_places, "f"), bootstrap_conf_low),
               ", ",
               sprintf(paste0("%.", decimal_places, "f"), bootstrap_conf_high),
               ")")
      } else {
        NA
      },

      # Get sample size (with safety check for missing column)
      sample_size_col = if("sample_size" %in% names(.)) {
        ifelse(!is.na(.data$sample_size), .data$sample_size, "")
      } else {
        "N/A"
      },

      # Create balance ratio (as fraction for easy interpretation)
      balance_ratio = if(all(c("n_balanced_vars", "n_total_vars") %in% names(.))) {
        ifelse(!is.na(.data$n_balanced_vars) & !is.na(.data$n_total_vars) & .data$n_total_vars > 0,
               paste0(.data$n_balanced_vars, "/", .data$n_total_vars),
               NA)
      } else {
        NA
      },

      # Create significance indicator
      is_significant = glm_p_value < alpha_threshold,

      # Create significance symbol
      sig_symbol = case_when(
        glm_p_value < 0.001 ~ "***",
        glm_p_value < 0.01 ~ "**",
        glm_p_value < 0.05 ~ "*",
        TRUE ~ ""
      ),

      # Format p-value with significance stars
      p_value_formatted = paste0(sprintf(paste0("%.", decimal_places, "f"), glm_p_value), sig_symbol)
    ) %>%
    dplyr::arrange(glm_p_value) %>%  # Sort by significance (most significant first)
    dplyr::mutate(
      row_id = row_number()
    )

  # Select and rename columns for the table
  if (include_bootstrap && "bootstrap_estimate" %in% names(formatted_data) && !all(is.na(formatted_data$bootstrap_estimate))) {
    table_final <- formatted_data %>%
      dplyr::select(
        row_id,
        pretty_model,
        glm_estimate,
        bootstrap_estimate,
        glm_conf_int,
        bootstrap_conf_int,
        p_value_formatted,
        balance_ratio,
        sample_size_col,
        is_significant
      )

    col_labels <- list(
      pretty_model = "Model",
      glm_estimate = "Estimate",
      bootstrap_estimate = "Bootstrap Est.",
      glm_conf_int = "95% CI",
      bootstrap_conf_int = "Bootstrap 95% CI",
      p_value_formatted = "P-Value",
      balance_ratio = "Balanced",
      sample_size_col = "N"
    )
  } else {
    table_final <- formatted_data %>%
      dplyr::select(
        row_id,
        pretty_model,
        glm_estimate,
        glm_conf_int,
        p_value_formatted,
        balance_ratio,
        sample_size_col,
        is_significant
      )

    col_labels <- list(
      pretty_model = "Model",
      glm_estimate = "Estimate",
      glm_conf_int = "95% CI",
      p_value_formatted = "P-Value",
      balance_ratio = "Balanced",
      sample_size_col = "N"
    )
  }

  # Create the GT table
  gt_table <- table_final %>%
    gt::gt() %>%
    gt::cols_label(.list = col_labels) %>%
    gt::fmt_number(
      columns = c("glm_estimate",
                  if("bootstrap_estimate" %in% names(table_final)) "bootstrap_estimate" else NULL),
      decimals = decimal_places
    ) %>%
    gt::tab_header(
      title = gt::md(table_title),
      subtitle = gt::md(paste("Treatment Variable:", treatment_label))
    ) %>%

    # Style the model names
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_body(columns = "pretty_model")
    ) %>%

    # Highlight significant results
    gt::tab_style(
      style = list(
        gt::cell_fill(color = "#e8f5e8"),
        gt::cell_text(weight = "bold")
      ),
      locations = gt::cells_body(rows = table_final$row_id[table_final$is_significant])
    ) %>%

    # Style non-significant results
    gt::tab_style(
      style = gt::cell_text(color = "#666666"),
      locations = gt::cells_body(rows = table_final$row_id[!table_final$is_significant])
    ) %>%

    # Set font size
    gt::tab_style(
      style = gt::cell_text(size = font_size),
      locations = list(
        gt::cells_column_labels(columns = gt::everything()),
        gt::cells_body(columns = gt::everything())
      )
    ) %>%

    # Table options
    gt::tab_options(
      table.font.size = gt::px(font_size),
      table.font.names = font_family,
      data_row.padding = gt::px(8),
      row_group.padding = gt::px(8),
      heading.padding = gt::px(8)
    ) %>%

    # Alignment
    gt::cols_align(align = "left", columns = "pretty_model") %>%
    gt::cols_align(align = "center", columns = c("sample_size_col", "balance_ratio", "p_value_formatted")) %>%
    gt::cols_align(align = "right", columns = contains("estimate")) %>%

    # Hide utility columns
    gt::cols_hide(columns = c("row_id", "is_significant")) %>%

    # Add footnote explaining significance stars
    gt::tab_footnote(
      footnote = "* p < 0.05, ** p < 0.01, *** p < 0.001",
      locations = gt::cells_column_labels(columns = "p_value_formatted")
    ) %>%

    # Add footnote about highlighting
    gt::tab_footnote(
      footnote = paste("Significant results (p <", alpha_threshold, ") are highlighted in green"),
      locations = gt::cells_title()
    )

  # Add balance footnote if balance data is available
  if ("balance_ratio" %in% names(table_final) && any(!is.na(table_final$balance_ratio))) {
    gt_table <- gt_table %>%
      gt::tab_footnote(
        footnote = "Balanced = Number of covariates with |standardized difference| < 0.1",
        locations = gt::cells_column_labels(columns = "balance_ratio")
      )
  }

  # Save if filename provided
  if (!is.null(filename)) {
    gt::gtsave(gt_table, filename = filename)
    cat("Table saved as:", filename, "\n")
  }

  # Print summary
  sig_count <- sum(table_final$is_significant)
  total_count <- nrow(table_final)
  cat("\nTable Summary:\n")
  cat("  Total models:", total_count, "\n")
  cat("  Significant models (p <", alpha_threshold, "):", sig_count, "\n")
  cat("  Treatment variable:", treatment_label, "\n")

  return(gt_table)
}


#' Quick Treatment Table
#'
#' Creates a quick table showing only significant treatment effects.
#'
#' @param sig_results Output from \code{\link{check_treatment_significance}}
#' @param title Character. Table title
#' @param treatment_name Character. Treatment variable name for labeling
#'
#' @return A gt table object
#'
#' @details
#' Convenience function that creates a table showing only models with significant
#' treatment effects. Equivalent to calling \code{create_significance_table}
#' with \code{show_all_models = FALSE}.
#'
#' @examples
#' \dontrun{
#' quick_treatment_table(sig_results, "Significant Treatment Effects")
#' }
#'
#' @seealso \code{\link{create_significance_table}}, \code{\link{full_treatment_table}}
#' @export
quick_treatment_table <- function(sig_results, title = "Treatment Significance Across Models", treatment_name = NULL) {
  create_significance_table(
    sig_results = sig_results,
    table_title = title,
    treatment_label = treatment_name,
    include_bootstrap = TRUE,
    decimal_places = 3,
    show_all_models = FALSE,  # Only show significant models by default
    alpha_threshold = 0.05
  )
}


#' Full Treatment Table
#'
#' Creates a table showing treatment effects from all models (significant and non-significant).
#'
#' @param sig_results Output from \code{\link{check_treatment_significance}}
#' @param title Character. Table title
#' @param treatment_name Character. Treatment variable name for labeling
#'
#' @return A gt table object
#'
#' @details
#' Convenience function that creates a table showing treatment effects from all models,
#' regardless of significance. Equivalent to calling \code{create_significance_table}
#' with \code{show_all_models = TRUE}.
#'
#' @examples
#' \dontrun{
#' full_treatment_table(sig_results, "All Treatment Effects")
#' }
#'
#' @seealso \code{\link{create_significance_table}}, \code{\link{quick_treatment_table}}
#' @export
full_treatment_table <- function(sig_results, title = "Treatment Effects: All Models", treatment_name = NULL) {
  create_significance_table(
    sig_results = sig_results,
    table_title = title,
    treatment_label = treatment_name,
    include_bootstrap = TRUE,
    decimal_places = 3,
    show_all_models = TRUE,  # Show all models
    alpha_threshold = 0.05
  )
}
