#' Create Publication-Ready Covariate Balance Table
#'
#' Creates a formatted table of covariate balance statistics across one or more
#' CBPS model specifications, suitable for publication or presentation.
#'
#' @param balance_results List or data frame. Either:
#'   \itemize{
#'     \item Single output from \code{\link{extract_cbps_results}} with include_balance_details = TRUE
#'     \item Named list of multiple outputs from \code{\link{extract_cbps_results}}
#'     \item Data frame with columns: variable, diff_adj, var_ratio_adj, and optionally model
#'   }
#' @param table_title Character. Title for the table (default: "Covariate Balance After Weighting")
#' @param subtitle Character. Optional subtitle for the table (default: NULL)
#' @param filename Character. Optional filename to save the table (default: NULL)
#' @param model_labels Named character vector. Pretty names for models (default: NULL)
#' @param variable_labels Named character vector or function. Pretty names for variables.
#'   If a function, it should take variable names and return formatted names (default: NULL)
#' @param decimal_places Integer. Number of decimal places for balance values (default: 3)
#' @param font_size Integer. Font size for the table (default: 16)
#' @param font_family Character. Font family for the table (default: "Times New Roman")
#' @param shade_alternating Logical. Whether to shade alternating model groups (default: TRUE)
#'
#' @return A gt table object that can be displayed or saved
#'
#' @details
#' This function creates publication-ready covariate balance tables using the gt package with:
#' \itemize{
#'   \item Professional formatting suitable for academic papers
#'   \item Both SMD (standardized mean difference) and VR (variance ratio) columns
#'   \item Alternating row shading by model groups for easy reading
#'   \item Support for single or multiple model comparisons
#'   \item Customizable variable and model labels
#' }
#'
#' @examples
#' \dontrun{
#' # Single model
#' results <- extract_cbps_results(cbps_output, include_balance_details = TRUE)
#' balance_table <- create_balance_table(
#'   balance_results = results,
#'   table_title = "Hypothesis 1: Parity Status",
#'   subtitle = "Covariate Balance After Weighting"
#' )
#'
#' # Multiple models
#' results_list <- list(
#'   "Model 1" = extract_cbps_results(cbps_output1, include_balance_details = TRUE),
#'   "Model 2" = extract_cbps_results(cbps_output2, include_balance_details = TRUE)
#' )
#' balance_table <- create_balance_table(
#'   balance_results = results_list,
#'   table_title = "Hypothesis 1: Parity Status",
#'   subtitle = "Covariate Balance After Weighting",
#'   model_labels = c("Model 1" = "Basic Specification", "Model 2" = "Full Model")
#' )
#' }
#'
#' @seealso \code{\link{extract_cbps_results}}
#' @export
#' @importFrom gt gt cols_label fmt_number tab_header tab_style cell_text cell_fill cells_body cells_column_labels cols_hide cols_align tab_options gtsave md px
#' @importFrom dplyr select mutate filter bind_rows row_number
#' @importFrom rlang .data
create_balance_table <- function(
    balance_results,
    table_title = "Covariate Balance After Weighting",
    subtitle = NULL,
    filename = NULL,
    model_labels = NULL,
    variable_labels = NULL,
    decimal_places = 3,
    font_size = 16,
    font_family = "Times New Roman",
    shade_alternating = TRUE,
    smd_threshold = 0.1,
    vr_lower = 0.5,
    vr_upper = 2.0,
    show_model_column = NULL  # Changed to NULL default - auto-detect
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

    # Check for balance columns (various possible names)
    if ("diff_adj" %in% names(balance_data)) {
      balance_data <- balance_data %>%
        dplyr::mutate(mean_balance = abs(.data$diff_adj))
    } else if ("mean_balance" %in% names(balance_data)) {
      # Already has mean_balance, keep it
      balance_data <- balance_data
    } else {
      stop("Input data frame must contain 'diff_adj' or 'mean_balance' column")
    }

    # Check for variance ratio
    if ("var_ratio_adj" %in% names(balance_data)) {
      balance_data <- balance_data %>%
        dplyr::mutate(variance_balance = .data$var_ratio_adj)
    } else if (!"variance_balance" %in% names(balance_data)) {
      # If no variance ratio exists, create NA column
      balance_data$variance_balance <- NA
    }

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
          mean_balance = abs(.data$diff_adj),
          variance_balance = if("var_ratio_adj" %in% names(.)) .data$var_ratio_adj else NA
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
              mean_balance = abs(.data$diff_adj),
              variance_balance = if("var_ratio_adj" %in% names(.)) .data$var_ratio_adj else NA
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

  # Check number of models
  n_models <- length(unique(balance_data$model))

  # STEP 2: Apply pretty labels for models
  if (!is.null(model_labels)) {
    balance_data <- balance_data %>%
      dplyr::mutate(
        model = ifelse(.data$model %in% names(model_labels),
                       model_labels[.data$model],
                       .data$model)
      )
  }

  # FIXED: Determine if model column should be hidden
  # Logic: Hide if (1) only one model AND (2) no custom model labels provided AND (3) show_model_column not explicitly set to TRUE
  if (is.null(show_model_column)) {
    # Auto-detect: hide if single model with no custom labels
    hide_model_column <- (n_models == 1 && is.null(model_labels))
  } else {
    # User explicitly set preference
    hide_model_column <- !show_model_column
  }

  # STEP 3: Apply pretty labels for variables
  if (!is.null(variable_labels)) {
    if (is.function(variable_labels)) {
      # Apply function to create pretty labels
      balance_data <- balance_data %>%
        dplyr::mutate(pretty_variable = variable_labels(.data$variable))
    } else if (is.character(variable_labels)) {
      # Use named vector for mapping
      balance_data <- balance_data %>%
        dplyr::mutate(
          pretty_variable = ifelse(.data$variable %in% names(variable_labels),
                                   variable_labels[.data$variable],
                                   .data$variable)
        )
    }
  } else {
    # No labels provided - use variable names as-is
    balance_data <- balance_data %>%
      dplyr::mutate(pretty_variable = .data$variable)
  }

  # STEP 3.5: Flag poor balance and add asterisk
  balance_data <- balance_data %>%
    dplyr::mutate(
      poor_balance = (.data$mean_balance > smd_threshold) |
        (.data$variance_balance < vr_lower) |
        (.data$variance_balance > vr_upper),
      # Add asterisk to variable name if poorly balanced
      pretty_variable = ifelse(.data$poor_balance,
                               paste0(.data$pretty_variable, "*"),
                               .data$pretty_variable)
    )

  # Count number of poorly balanced covariates
  n_poor_balance <- sum(balance_data$poor_balance, na.rm = TRUE)

  # STEP 4: Format the data for table
  table_data <- balance_data %>%
    dplyr::mutate(
      model = as.character(.data$model),
      model_display = ifelse(duplicated(.data$model), "", .data$model)
    ) %>%
    dplyr::filter(!is.na(.data$variance_balance) | !is.na(.data$mean_balance)) %>%
    # Assign a block ID to each group of rows by model
    dplyr::mutate(group_id = with(rle(.data$model), rep(seq_along(values), lengths))) %>%
    # Create shading flag
    dplyr::mutate(
      shade_flag = shade_alternating & (.data$group_id %% 2 == 1),
      row_id = dplyr::row_number()
    ) %>%
    dplyr::select(.data$row_id, .data$model_display, .data$pretty_variable,
                  .data$mean_balance, .data$variance_balance, .data$shade_flag, .data$poor_balance)

  # STEP 5: Create the gt table
  gt_table <- table_data %>%
    gt::gt() %>%
    gt::cols_label(
      model_display = "Model",
      pretty_variable = "Covariate",
      mean_balance = "SMD",
      variance_balance = "VR"
    ) %>%
    gt::fmt_number(
      columns = c("mean_balance", "variance_balance"),
      decimals = decimal_places
    )

  # Add header with optional subtitle
  if (!is.null(subtitle)) {
    gt_table <- gt_table %>%
      gt::tab_header(
        title = gt::md(table_title),
        subtitle = gt::md(subtitle)
      )
  } else {
    gt_table <- gt_table %>%
      gt::tab_header(
        title = gt::md(table_title)
      )
  }

  # Apply styling
  gt_table <- gt_table %>%
    # Bold the model display values (non-empty cells)
    gt::tab_style(
      style = gt::cell_text(weight = "bold"),
      locations = gt::cells_body(
        columns = "model_display",
        rows = table_data$model_display != ""
      )
    ) %>%

    # Set font size for all elements
    gt::tab_style(
      style = gt::cell_text(size = font_size),
      locations = list(
        gt::cells_column_labels(columns = gt::everything()),
        gt::cells_body(columns = gt::everything())
      )
    ) %>%

    # Apply highlight to every other model block
    gt::tab_style(
      style = gt::cell_fill(color = "#f7f7f7"),
      locations = gt::cells_body(
        rows = table_data$row_id[table_data$shade_flag]
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
    gt::cols_align(align = "left", columns = c("pretty_variable", "model_display")) %>%
    gt::cols_align(align = "center", columns = c("mean_balance", "variance_balance"))

  # Hide utility columns and optionally model column
  columns_to_hide <- c("row_id", "shade_flag", "poor_balance")
  if (hide_model_column) {
    columns_to_hide <- c(columns_to_hide, "model_display")
  }

  gt_table <- gt_table %>%
    gt::cols_hide(columns = columns_to_hide) %>%

    # Add source note explaining asterisk
    gt::tab_source_note(
      source_note = gt::md(paste0(
        "*Indicates poor balance (SMD > ", smd_threshold,
        " or VR < ", vr_lower, " or VR > ", vr_upper, ")"
      ))
    )

  # Save if filename provided
  if (!is.null(filename)) {
    gt::gtsave(gt_table, filename = filename)
    cat("Table saved as:", filename, "\n")
  }

  # Print summary
  n_total <- nrow(table_data)

  cat("\nBalance Table Summary:\n")
  cat("  Number of models:", n_models, "\n")
  cat("  Total covariate rows:", n_total, "\n")
  cat("  Poorly balanced covariates:", n_poor_balance, "\n")
  if (hide_model_column) {
    cat("  Model column hidden (single model, no custom labels)\n")
  }

  return(gt_table)
}
