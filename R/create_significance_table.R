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
#' @param include_fit_stats Logical. Whether to include model fit statistics (AIC/BIC/BICc) (default: TRUE)
#' @param decimal_places Integer. Number of decimal places for estimates (default: 3)
#' @param font_size Integer. Font size for the table (default: 16)
#' @param font_family Character. Font family for the table (default: "Times New Roman")
#' @param show_all_models Logical. Whether to show all models or just significant ones (default: FALSE)
#' @param alpha_threshold Numeric. Threshold for highlighting significance (default: 0.05)
#' @param include_subtitle Logical. Whether to include treatment variable subtitle (default: FALSE)
#'
#' @return A gt table object that can be displayed or saved
#'
#' @details
#' This function creates publication-ready tables using the gt package with:
#' \itemize{
#'   \item Professional formatting suitable for academic papers
#'   \item Significance highlighting (light blue background for significant results)
#'   \item Confidence intervals for both GLM and bootstrap estimates
#'   \item Balance assessment ratios
#'   \item Model fit statistics (AIC, BIC, BICc)
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
#' @seealso \code{\link{check_treatment_significance}}
#' @export
#' @importFrom gt gt cols_label fmt_number tab_header tab_style cell_text cell_fill cells_body cells_column_labels cells_title cols_hide cols_align tab_footnote tab_options gtsave md px
#' @importFrom dplyr select all_of mutate case_when arrange filter
create_significance_table <- function(
    sig_results,
    table_title = "Treatment Effects Across Model Specifications",
    filename = NULL,
    treatment_label = NULL,
    model_labels = NULL,
    include_fit_stats = TRUE,
    include_mi_info = TRUE,
    decimal_places = 3,
    font_size = 16,
    font_family = "Times New Roman",
    show_all_models = FALSE,
    alpha_threshold = 0.05,
    include_subtitle = FALSE
) {

  # Check required packages
  required_packages <- c("dplyr", "gt", "tibble")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

  if (length(missing_packages) > 0) {
    stop("Missing required packages: ", paste(missing_packages, collapse = ", "),
         "\nTry restarting R and reinstalling CBPSAnalysis package.")
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
    treatment_label <- unique(table_data$variable)[1]
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

      # Format confidence intervals
      glm_conf_int = if(all(c("glm_conf_low", "glm_conf_high") %in% names(.))) {
        paste0("(",
               sprintf(paste0("%.", decimal_places, "f"), glm_conf_low),
               ", ",
               sprintf(paste0("%.", decimal_places, "f"), glm_conf_high),
               ")")
      } else {
        "CI not available"
      },

      # Get sample size - handle both sample_size and final_sample_size columns
      sample_size_col = if("sample_size" %in% names(.)) {
        ifelse(!is.na(.data$sample_size), .data$sample_size, "")
      } else if("final_sample_size" %in% names(.)) {
        ifelse(!is.na(.data$final_sample_size), .data$final_sample_size, "")
      } else {
        "N/A"
      },

      # Format MI info if available
      mi_info = if(include_mi_info && "n_imputations" %in% names(.) && "inference_method" %in% names(.)) {
        ifelse(.data$n_imputations > 1,
               paste0("MI (m=", .data$n_imputations, ")"),
               "Single")
      } else {
        NA
      },

      # Format FMI if available
      fmi_formatted = if(include_mi_info && "fmi" %in% names(.)) {
        ifelse(!is.na(.data$fmi),
               sprintf(paste0("%.", decimal_places, "f"), .data$fmi),
               NA)
      } else {
        NA
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
    dplyr::arrange(glm_p_value) %>%
    dplyr::mutate(
      row_id = row_number()
    )

  # Select and rename columns for the table (NO BOOTSTRAP)
  base_select <- c(
    "row_id",
    "pretty_model",
    "glm_estimate",
    "glm_conf_int",
    "p_value_formatted",
    "balance_ratio",
    "sample_size_col",
    "is_significant"
  )

  # Add MI info columns if they exist and aren't all NA
  if (include_mi_info && "mi_info" %in% names(formatted_data) && !all(is.na(formatted_data$mi_info))) {
    base_select <- c(base_select, "mi_info")
  }
  if (include_mi_info && "fmi_formatted" %in% names(formatted_data) && !all(is.na(formatted_data$fmi_formatted))) {
    base_select <- c(base_select, "fmi_formatted")
  }

  # Add fit stats if requested and available
  if (include_fit_stats && all(c("AIC", "BIC", "BICc") %in% names(formatted_data))) {
    base_select <- c(base_select, "AIC", "BIC", "BICc")
  }

  table_final <- formatted_data %>%
    dplyr::select(dplyr::all_of(base_select))

  col_labels <- list(
    pretty_model = "Model",
    glm_estimate = "Estimate",
    glm_conf_int = "95% CI",
    p_value_formatted = "P-Value",
    balance_ratio = "Covariates Balanced",
    sample_size_col = "N"
  )

  # Add MI column labels if present
  if ("mi_info" %in% names(table_final)) {
    col_labels$mi_info <- "Method"
  }
  if ("fmi_formatted" %in% names(table_final)) {
    col_labels$fmi_formatted <- "FMI"
  }

  if (include_fit_stats && "AIC" %in% names(table_final)) {
    col_labels$AIC <- "AIC"
    col_labels$BIC <- "BIC"
    col_labels$BICc <- "BICc"
  }

  # Create the GT table with conditional subtitle
  if (include_subtitle) {
    gt_table <- table_final %>%
      gt::gt() %>%
      gt::cols_label(.list = col_labels) %>%
      gt::tab_header(
        title = gt::md(table_title),
        subtitle = gt::md(paste("Treatment Variable:", treatment_label))
      )
  } else {
    gt_table <- table_final %>%
      gt::gt() %>%
      gt::cols_label(.list = col_labels) %>%
      gt::tab_header(
        title = gt::md(table_title)
      )
  }

  # Format numeric columns
  gt_table <- gt_table %>%
    gt::fmt_number(
      columns = "glm_estimate",
      decimals = decimal_places
    )

  # Format fit statistics with 1 decimal place if present
  if (include_fit_stats && all(c("AIC", "BIC", "BICc") %in% names(table_final))) {
    gt_table <- gt_table %>%
      gt::fmt_number(
        columns = c("AIC", "BIC", "BICc"),
        decimals = 1
      )
  }

  gt_table <- gt_table %>%
    # Style the model names
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_body(columns = "pretty_model")
    ) %>%

    # Highlight significant results with light blue
    gt::tab_style(
      style = list(
        gt::cell_fill(color = "#e3f2fd"),
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
    gt::cols_align(align = "right", columns = "glm_estimate")

  # Align MI columns to center if present
  if ("mi_info" %in% names(table_final)) {
    gt_table <- gt_table %>%
      gt::cols_align(align = "center", columns = "mi_info")
  }
  if ("fmi_formatted" %in% names(table_final)) {
    gt_table <- gt_table %>%
      gt::cols_align(align = "center", columns = "fmi_formatted")
  }

  # Align fit statistics to center if present
  if (include_fit_stats && "AIC" %in% names(table_final)) {
    gt_table <- gt_table %>%
      gt::cols_align(align = "center", columns = c("AIC", "BIC", "BICc"))
  }

  gt_table <- gt_table %>%
    # Hide utility columns
    gt::cols_hide(columns = c("row_id", "is_significant"))

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

  # Show MI info in summary
  if ("mi_info" %in% names(table_final)) {
    mi_models <- sum(grepl("MI", table_final$mi_info, ignore.case = TRUE))
    if (mi_models > 0) {
      cat("  Models using multiple imputation:", mi_models, "\n")
    }
  }

  if (include_fit_stats && "AIC" %in% names(table_final)) {
    best_aic <- table_final %>% dplyr::slice_min(AIC, n = 1)
    cat("  Best fitting model (lowest AIC):", best_aic$pretty_model[1], "\n")
  }

  return(gt_table)
}
