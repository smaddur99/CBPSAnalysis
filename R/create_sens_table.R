#' Create Sensitivity Analysis Results Table
#'
#' Creates a formatted table displaying treatment effects across multiple sensitivity analyses.
#' The table is organized with sensitivity analysis types as row groups, allowing comparison
#' of results across different model specifications (e.g., main model, adjusted models).
#'
#' @param sig_results Data frame or list containing sensitivity analysis results.
#'   Must include columns: variable, model_name, glm_estimate, glm_p_value, sensitivity_type.
#'   Optional columns: glm_conf_low, glm_conf_high, bootstrap_estimate, bootstrap_conf_low,
#'   bootstrap_conf_high, sample_size.
#'   Use \code{prepare_sens_table_data()} to properly format output from multiple
#'   \code{cbps_weighted_analysis()} runs.
#' @param table_title Character. Title for the table (default: "Sensitivity Analysis Results")
#' @param filename Character or NULL. Filename to save table (e.g., "sens_table.html").
#'   If NULL, table is not saved (default: NULL)
#' @param treatment_label Character or NULL. Label for treatment variable. If NULL, uses first
#'   unique variable name from data (default: NULL)
#' @param sens_labels Named character vector or NULL. Custom labels for sensitivity analysis types.
#'   Names should match sensitivity_type values in sig_results. If NULL, uses original
#'   sensitivity_type values (default: NULL)
#' @param include_bootstrap Logical. Whether to include bootstrap estimates and confidence
#'   intervals in table (default: TRUE)
#' @param decimal_places Integer. Number of decimal places for estimates (default: 3)
#' @param font_size Integer. Font size in pixels (default: 16)
#' @param font_family Character. Font family for table (default: "Times New Roman")
#' @param show_all_models Logical. If sig_results contains both 'significant_only' and
#'   'detailed_results', whether to show all models (TRUE) or only significant ones (FALSE).
#'   (default: FALSE)
#' @param alpha_threshold Numeric. Significance threshold for highlighting rows (default: 0.05)
#' @param include_subtitle Logical. Whether to include treatment variable in subtitle (default: FALSE)
#'
#' @return A gt table object with formatted sensitivity analysis results
#'
#' @details
#' The function creates a publication-ready table organized by sensitivity analysis type.
#' Each sensitivity type appears as a row group, with outcomes listed within each group.
#' Significant results (p < alpha_threshold) are highlighted with blue background and bold text.
#' Non-significant results appear in gray.
#'
#' The table includes:
#' \itemize{
#'   \item Point estimates (from GLM and optionally bootstrap)
#'   \item 95% confidence intervals
#'   \item P-values with significance stars (*** p<0.001, ** p<0.01, * p<0.05)
#'   \item Sample sizes
#' }
#'
#' Input data should typically be prepared using \code{prepare_sens_table_data()}, which
#' extracts and formats results from multiple \code{cbps_weighted_analysis()} outputs.
#'
#' @seealso \code{\link{prepare_sens_table_data}} for preparing input data,
#'   \code{\link{create_f_sig_table}} for single model specification tables
#'
#' @examples
#' \dontrun{
#' # Prepare data from multiple model runs
#' sens_data <- prepare_sens_table_data(
#'   model_list = list(
#'     main = list(outcome1 = main_model1, outcome2 = main_model2),
#'     adj_z = list(outcome1 = sens_model1, outcome2 = sens_model2)
#'   ),
#'   sens_labels = c(main = "Main Model", adj_z = "Adjusted for Z"),
#'   treatment_var = "pregnancy_advice"
#' )
#'
#' # Create sensitivity analysis table
#' create_sens_table(
#'   sens_data,
#'   table_title = "Sensitivity Analysis: Treatment Effects",
#'   treatment_label = "Pregnancy Advice",
#'   include_bootstrap = TRUE,
#'   filename = "sensitivity_table.html"
#' )
#'
#' # Customize appearance
#' create_sens_table(
#'   sens_data,
#'   sens_labels = c(
#'     main = "Main Model (Unadjusted)",
#'     adj_z = "Adjusted for Maternal Age",
#'     adj_z_income = "Adjusted for Age + Income"
#'   ),
#'   decimal_places = 2,
#'   font_size = 14,
#'   alpha_threshold = 0.01
#' )
#' }
#'
#' @export
#' @importFrom dplyr filter mutate select arrange group_by case_when
#' @importFrom gt gt cols_label tab_header fmt_number tab_style cell_text cell_fill cells_body cells_column_labels tab_options cols_align cols_hide md
#' @importFrom tibble tibble
create_sens_table <- function(sig_results,
                              table_title = "Sensitivity Analysis Results",
                              filename = NULL,
                              treatment_label = NULL,
                              sens_labels = NULL,
                              include_bootstrap = TRUE,
                              decimal_places = 3,
                              font_size = 16,
                              font_family = "Times New Roman",
                              show_all_models = FALSE,
                              alpha_threshold = 0.05,
                              include_subtitle = FALSE) {

  required_packages <- c("dplyr", "gt", "tibble")
  missing_packages <- required_packages[!sapply(required_packages,
                                                requireNamespace, quietly = TRUE)]
  if (length(missing_packages) > 0) {
    stop("Missing required packages: ", paste(missing_packages, collapse = ", "))
  }

  # Extract data
  if ("significant_only" %in% names(sig_results)) {
    if (show_all_models) {
      table_data <- sig_results$detailed_results
    } else {
      table_data <- sig_results$significant_only
    }
  } else {
    table_data <- sig_results
  }

  # Check for required columns
  required_cols <- c("variable", "glm_estimate", "glm_p_value",
                     "model_name", "sensitivity_type")
  missing_cols <- required_cols[!required_cols %in% names(table_data)]
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "),
         "\nDid you use prepare_sens_table_data() to format your results?")
  }

  # Set default labels
  if (is.null(treatment_label)) {
    treatment_label <- unique(table_data$variable)[1]
  }

  if (is.null(sens_labels)) {
    sens_labels <- setNames(unique(table_data$sensitivity_type),
                            unique(table_data$sensitivity_type))
  }

  # Format data
  formatted_data <- table_data %>%
    dplyr::mutate(
      sens_group = ifelse(sensitivity_type %in% names(sens_labels),
                          sens_labels[sensitivity_type],
                          sensitivity_type),
      outcome_name = model_name,
      glm_conf_int = if (all(c("glm_conf_low", "glm_conf_high") %in% names(.))) {
        paste0("(", sprintf(paste0("%.", decimal_places, "f"), glm_conf_low),
               ", ", sprintf(paste0("%.", decimal_places, "f"), glm_conf_high), ")")
      } else {
        "CI not available"
      },
      bootstrap_conf_int = if (include_bootstrap &&
                               all(c("bootstrap_conf_low", "bootstrap_conf_high") %in% names(.))) {
        paste0("(", sprintf(paste0("%.", decimal_places, "f"), bootstrap_conf_low),
               ", ", sprintf(paste0("%.", decimal_places, "f"), bootstrap_conf_high), ")")
      } else {
        NA
      },
      # Get sample size - handle both sample_size and final_sample_size columns
      sample_size_col = if("sample_size" %in% names(.)) {
        ifelse(!is.na(.data$sample_size), .data$sample_size, "")
      } else if("final_sample_size" %in% names(.)) {
        ifelse(!is.na(.data$final_sample_size), .data$final_sample_size, "")
      } else {
        "N/A"
      },
      is_significant = glm_p_value < alpha_threshold,
      sig_symbol = case_when(
        glm_p_value < 0.001 ~ "***",
        glm_p_value < 0.01 ~ "**",
        glm_p_value < 0.05 ~ "*",
        TRUE ~ ""
      ),
      p_value_formatted = paste0(sprintf(paste0("%.", decimal_places, "f"),
                                         glm_p_value), sig_symbol)
    ) %>%
    dplyr::arrange(sensitivity_type, glm_p_value) %>%
    dplyr::mutate(row_id = row_number())

  # Select columns - KEEP sens_group in data but DON'T label it
  if (include_bootstrap && "bootstrap_estimate" %in% names(formatted_data) &&
      !all(is.na(formatted_data$bootstrap_estimate))) {
    table_final <- formatted_data %>%
      dplyr::select(row_id, sens_group, outcome_name, glm_estimate,
                    bootstrap_estimate, glm_conf_int, bootstrap_conf_int,
                    p_value_formatted, sample_size_col, is_significant)

    # CRITICAL: Do NOT include sens_group in col_labels
    col_labels <- list(
      outcome_name = "Outcome",
      glm_estimate = "Estimate",
      bootstrap_estimate = "Bootstrap Est.",
      glm_conf_int = "95% CI",
      bootstrap_conf_int = "Bootstrap 95% CI",
      p_value_formatted = "P-Value",
      sample_size_col = "N"
    )
  } else {
    table_final <- formatted_data %>%
      dplyr::select(row_id, sens_group, outcome_name, glm_estimate,
                    glm_conf_int, p_value_formatted, sample_size_col, is_significant)

    # CRITICAL: Do NOT include sens_group in col_labels
    col_labels <- list(
      outcome_name = "Outcome",
      glm_estimate = "Estimate",
      glm_conf_int = "95% CI",
      p_value_formatted = "P-Value",
      sample_size_col = "N"
    )
  }

  # Create gt table with row groups
  if (include_subtitle) {
    gt_table <- table_final %>%
      gt::gt(groupname_col = "sens_group") %>%  # CHANGED: Use groupname_col parameter
      gt::cols_label(.list = col_labels) %>%
      gt::tab_header(
        title = gt::md(table_title),
        subtitle = gt::md(paste("Treatment Variable:", treatment_label))
      )
  } else {
    gt_table <- table_final %>%
      gt::gt(groupname_col = "sens_group") %>%  # CHANGED: Use groupname_col parameter
      gt::cols_label(.list = col_labels) %>%
      gt::tab_header(title = gt::md(table_title))
  }

  # Apply styling
  gt_table <- gt_table %>%
    gt::fmt_number(
      columns = c("glm_estimate",
                  if ("bootstrap_estimate" %in% names(table_final)) "bootstrap_estimate" else NULL),
      decimals = decimal_places
    ) %>%
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_body(columns = "outcome_name")
    ) %>%
    gt::tab_style(
      style = list(gt::cell_fill(color = "#e3f2fd"),
                   gt::cell_text(weight = "bold")),
      locations = gt::cells_body(rows = table_final$row_id[table_final$is_significant])
    ) %>%
    gt::tab_style(
      style = gt::cell_text(color = "#666666"),
      locations = gt::cells_body(rows = table_final$row_id[!table_final$is_significant])
    ) %>%
    gt::tab_style(
      style = gt::cell_text(size = font_size),
      locations = list(
        gt::cells_column_labels(columns = gt::everything()),
        gt::cells_body(columns = gt::everything())
      )
    ) %>%
    gt::tab_options(
      table.font.size = gt::px(font_size),
      table.font.names = font_family,
      data_row.padding = gt::px(8),
      row_group.padding = gt::px(8),
      heading.padding = gt::px(8)
    ) %>%
    gt::cols_align(align = "left", columns = "outcome_name") %>%
    gt::cols_align(align = "center", columns = c("sample_size_col", "p_value_formatted")) %>%
    gt::cols_align(align = "right", columns = contains("estimate")) %>%
    gt::cols_hide(columns = c("row_id", "is_significant"))  # Hide helper columns but NOT sens_group

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
  cat("  Number of sensitivity analyses:", length(unique(table_final$sens_group)), "\n")

  return(gt_table)
}
