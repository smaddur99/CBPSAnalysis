#' Create Publication-Ready Covariate Balance Table
#'
#' Creates a formatted table of covariate balance statistics across one or more
#' CBPS model specifications, suitable for publication or presentation.
#'
#' @param balance_results List or data frame. Either:
#'   \itemize{
#'     \item Single output from \code{\link{extract_cbps_results}} with include_balance_details = TRUE
#'     \item Named list of multiple outputs from \code{\link{extract_cbps_results}}
#'     \item Data frame with columns: variable, diff_adj (or mean_balance), and optionally model
#'   }
#' @param table_title Character. Title for the table (default: "Covariate Balance After Weighting")
#' @param subtitle Character. Optional subtitle for the table (default: NULL)
#' @param filename Character. Optional filename to save the table (default: NULL)
#' @param model_labels Named character vector. Pretty names for models (default: NULL)
#' @param variable_labels Named character vector or function. Pretty names for variables.
#'   If a function, it should take variable names and return formatted names (default: NULL)
#' @param smd_threshold Numeric. Threshold for flagging imbalanced covariates (default: 0.1)
#' @param decimal_places Integer. Number of decimal places for SMD values (default: 3)
#' @param font_size Integer. Font size for the table (default: 16)
#' @param font_family Character. Font family for the table (default: "Times New Roman")
#' @param show_threshold_flag Logical. Whether to add asterisk to values above threshold (default: TRUE)
#' @param shade_alternating Logical. Whether to shade alternating model groups (default: TRUE)
#'
#' @return A gt table object that can be displayed or saved
#'
#' @details
#' This function creates publication-ready covariate balance tables using the gt package with:
#' \itemize{
#'   \item Professional formatting suitable for academic papers
#'   \item Alternating row shading by model groups for easy reading
#'   \item Automatic flagging of imbalanced covariates (SMD > threshold)
#'   \item Support for single or multiple model comparisons
#'   \item Customizable variable and model labels
#' }
#'
#' The function works with output from \code{extract_cbps_results()} and automatically
#' extracts the balance statistics. It can handle:
#' \itemize{
#'   \item Single model results
#'   \item Multiple models for comparison
#'   \item Custom variable naming functions
#'   \item Different SMD thresholds
#' }
#'
#' @examples
#' \dontrun{
#' # Single model
#' results <- extract_cbps_results(cbps_output, include_balance_details = TRUE)
#' balance_table <- create_balance_table(
#'   balance_results = results,
#'   table_title = "Covariate Balance After CBPS Weighting"
#' )
#'
#' # Multiple models
#' results_list <- list(
#'   "Model 1" = extract_cbps_results(cbps_output1, include_balance_details = TRUE),
#'   "Model 2" = extract_cbps_results(cbps_output2, include_balance_details = TRUE)
#' )
#' balance_table <- create_balance_table(
#'   balance_results = results_list,
#'   table_title = "Comparison of Covariate Balance",
#'   model_labels = c("Model 1" = "Basic Specification", "Model 2" = "Full Model")
#' )
#'
#' # With custom variable labels
#' var_labels <- c(
#'   "age" = "Age (years)",
#'   "income" = "Annual Income ($)",
#'   "education" = "Education Level"
#' )
#' balance_table <- create_balance_table(
#'   balance_results = results,
#'   variable_labels = var_labels
#' )
#' }
#'
#' @seealso \code{\link{extract_cbps_results}}
#' @export
#' @importFrom gt gt cols_label fmt_number tab_header tab_style cell_text cell_fill cells_body cells_column_labels cols_hide cols_align tab_options gtsave md px
#' @importFrom dplyr select mutate case_when filter bind_rows arrange row_number
create_balance_table <- function(
    balance_results,
    table_title = "Covariate Balance After Weighting",
    subtitle = NULL,
    filename = NULL,
    model_labels = NULL,
    variable_labels = NULL,
    smd_threshold = 0.1,
    decimal_places = 3,
    font_size = 16,
    font_family = "Times New Roman",
    show_threshold_flag = TRUE,
    shade_alternating = TRUE
) {

  # Check required packages
  required_packages <- c("dplyr", "gt")
  missing_packages <- required_packages[!sapply(required_packages, requireNamespace, quietly = TRUE)]

  if (length(missing_packages) > 0) {
    stop("Missing required packages: ", paste(missing_packages, collapse = ", "),
         "\nPlease install these packages.")
  }

  # STEP 1: Extract and prepare balance data
  if (is.data.frame(balance_results)) {
    # Direct data frame input
    balance_data <- balance_results

    # Ensure required columns exist
    if (!"variable" %in% names(balance_data)) {
      stop("Input data frame must contain 'variable' column")
    }

    # Check for SMD column (various possible names)
    smd_col <- NULL
    if ("diff_adj" %in% names(balance_data)) {
      smd_col <- "diff_adj"
    } else if ("mean_balance" %in% names(balance_data)) {
      smd_col <- "mean_balance"
    } else if ("Diff.Adj" %in% names(balance_data)) {
      smd_col <- "Diff.Adj"
    } else {
      stop("Input data frame must contain a standardized mean difference column (diff_adj, mean_balance, or Diff.Adj)")
    }

    balance_data <- balance_data %>%
      dplyr::mutate(mean_balance = abs(!!rlang::sym(smd_col)))

    # Add model column if not present
    if (!"model" %in% names(balance_data)) {
      balance_data$model <- "Model"
    }

  } else if (is.list(balance_results)) {
    # Check if it's a single extract_cbps_results output
    if ("balance" %in% names(balance_results) && "estimates" %in% names(balance_results)) {
      # Single model from extract_cbps_results
      balance_data <- balance_results$balance %>%
        dplyr::mutate(
          model = "Model",
          mean_balance = abs(diff_adj)
        )
    } else {
      # Multiple models - named list of extract_cbps_results outputs
      balance_list <- list()

      for (model_name in names(balance_results)) {
        model_data <- balance_results[[model_name]]

        if ("balance" %in% names(model_data)) {
          temp_df <- model_data$balance %>%
            dplyr::mutate(
              model = model_name,
              mean_balance = abs(diff_adj)
            )
          balance_list[[model_name]] <- temp_df
        }
      }

      if (length(balance_list) == 0) {
        stop("No valid balance data found in the input list")
      }

      balance_data <- dplyr::bind_rows(balance_list)
    }
  } else {
    stop("balance_results must be a data frame or list output from extract_cbps_results()")
  }

  # STEP 2: Apply pretty labels for models
  if (!is.null(model_labels)) {
    balance_data <- balance_data %>%
      dplyr::mutate(
        model = ifelse(model %in% names(model_labels),
                       model_labels[model],
                       model)
      )
  }

  # STEP 3: Apply pretty labels for variables
  if (!is.null(variable_labels)) {
    if (is.function(variable_labels)) {
      # Apply function to create pretty labels
      balance_data <- balance_data %>%
        dplyr::mutate(pretty_variable = variable_labels(variable))
    } else if (is.character(variable_labels)) {
      # Use named vector for mapping
      balance_data <- balance_data %>%
        dplyr::mutate(
          pretty_variable = ifelse(variable %in% names(variable_labels),
                                   variable_labels[variable],
                                   variable)
        )
    }
  } else {
    # No labels provided - use variable names as-is
    balance_data <- balance_data %>%
      dplyr::mutate(pretty_variable = variable)
  }

  # STEP 4: Format the data for table
  table_data <- balance_data %>%
    dplyr::select(model, pretty_variable, mean_balance) %>%
    dplyr::filter(!is.na(mean_balance))

  # Add threshold flag if requested
  if (show_threshold_flag) {
    table_data <- table_data %>%
      dplyr::mutate(
        mean_balance_display = ifelse(
          mean_balance > smd_threshold,
          paste0(sprintf(paste0("%.", decimal_places, "f"), mean_balance), "*"),
          sprintf(paste0("%.", decimal_places, "f"), mean_balance)
        )
      )
  } else {
    table_data <- table_data %>%
      dplyr::mutate(
        mean_balance_display = sprintf(paste0("%.", decimal_places, "f"), mean_balance)
      )
  }

  # STEP 5: Prepare for alternating shading
  table_data <- table_data %>%
    dplyr::mutate(
      model_display = ifelse(duplicated(model), "", as.character(model))
    ) %>%
    dplyr::mutate(
      group_id = with(rle(model), rep(seq_along(values), lengths))
    ) %>%
    dplyr::mutate(
      shade_flag = shade_alternating & (group_id %% 2 == 1),
      row_id = dplyr::row_number()
    )

  # STEP 6: Create the gt table
  gt_table <- table_data %>%
    dplyr::select(row_id, model_display, pretty_variable, mean_balance_display, shade_flag) %>%
    gt::gt() %>%
    gt::cols_label(
      model_display = "Model",
      pretty_variable = "Covariate",
      mean_balance_display = "SMD"
    )

  # Add header with optional subtitle
  if (!is.null(subtitle)) {
    gt_table <- gt_table %>%
      gt::tab_header(
        title = gt::md(paste0("**", table_title, "**")),
        subtitle = gt::md(subtitle)
      )
  } else {
    gt_table <- gt_table %>%
      gt::tab_header(
        title = gt::md(paste0("**", table_title, "**"))
      )
  }

  # Apply styling
  gt_table <- gt_table %>%
    # Bold the model display values (non-empty cells)
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_body(
        columns = "model_display",
        rows = model_display != ""
      )
    ) %>%

    # General font styling
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
      heading.padding = gt::px(8)
    ) %>%

    # Alignment
    gt::cols_align(align = "left", columns = c("pretty_variable", "model_display")) %>%
    gt::cols_align(align = "center", columns = "mean_balance_display") %>%

    # Hide utility columns
    gt::cols_hide(columns = c("row_id", "shade_flag"))

  # Apply shading to alternating model groups if requested
  if (shade_alternating) {
    shaded_rows <- table_data$row_id[table_data$shade_flag]
    if (length(shaded_rows) > 0) {
      gt_table <- gt_table %>%
        gt::tab_style(
          style = gt::cell_fill(color = "#f7f7f7"),
          locations = gt::cells_body(rows = shaded_rows)
        )
    }
  }

  # Save if filename provided
  if (!is.null(filename)) {
    gt::gtsave(gt_table, filename = filename)
    cat("Table saved as:", filename, "\n")
  }

  # Print summary
  n_imbalanced <- sum(table_data$mean_balance > smd_threshold, na.rm = TRUE)
  n_total <- nrow(table_data)
  n_models <- length(unique(table_data$model))

  cat("\nBalance Table Summary:\n")
  cat("  Number of models:", n_models, "\n")
  cat("  Total covariates:", n_total, "\n")
  cat("  Imbalanced covariates (SMD >", smd_threshold, "):", n_imbalanced, "\n")
  if (show_threshold_flag) {
    cat("  Note: * indicates SMD >", smd_threshold, "\n")
  }

  return(gt_table)
}
