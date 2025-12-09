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
create_balance_table_grouped <- function(
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
    vr_upper = 2,
    use_row_groups = TRUE  # NEW: control whether to use row groups or column
) {

  # [Keep all your existing data processing code up to the table_data creation]
  # ... [all the validation and processing code stays the same] ...

  # Determine if we should use row groups
  n_models <- length(unique(balance_data$model))
  use_groups <- use_row_groups && (n_models > 1)

  if (use_groups) {
    # Prepare data for row groups (like sens_table)
    table_data <- balance_data %>%
      dplyr::mutate(
        model_group = ifelse(.data$model %in% names(model_labels),
                             model_labels[.data$model],
                             .data$model)
      ) %>%
      dplyr::filter(!is.na(.data$variance_balance) | !is.na(.data$mean_balance)) %>%
      dplyr::mutate(row_id = dplyr::row_number()) %>%
      dplyr::select(
        .data$row_id,
        .data$model_group,
        .data$pretty_variable,
        .data$mean_balance,
        .data$variance_balance,
        .data$poor_balance
      )

    # Create gt table with row groups
    gt_table <- table_data %>%
      gt::gt(groupname_col = "model_group") %>%
      gt::cols_label(
        pretty_variable = "Covariate",
        mean_balance = "SMD",
        variance_balance = "VR"
      ) %>%
      gt::fmt_number(
        columns = c("mean_balance", "variance_balance"),
        decimals = decimal_places
      )

  } else {
    # Original column-based approach
    table_data <- balance_data %>%
      dplyr::mutate(
        model = as.character(.data$model),
        model_display = ifelse(duplicated(.data$model), "", .data$model)
      ) %>%
      dplyr::filter(!is.na(.data$variance_balance) | !is.na(.data$mean_balance)) %>%
      dplyr::mutate(
        group_id = with(rle(.data$model), rep(seq_along(values), lengths))
      ) %>%
      dplyr::mutate(
        shade_flag = shade_alternating & (.data$group_id %% 2 == 1),
        row_id = dplyr::row_number()
      ) %>%
      dplyr::select(
        .data$row_id,
        .data$model_display,
        .data$pretty_variable,
        .data$mean_balance,
        .data$variance_balance,
        .data$shade_flag,
        .data$poor_balance
      )

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
      ) %>%
      gt::tab_style(
        style = gt::cell_fill(color = "#f7f7f7"),
        locations = gt::cells_body(rows = table_data$row_id[table_data$shade_flag])
      )
  }

  # [Keep all your existing styling code - it applies to both approaches]
  if (!is.null(subtitle)) {
    gt_table <- gt_table %>%
      gt::tab_header(title = gt::md(table_title), subtitle = gt::md(subtitle))
  } else {
    gt_table <- gt_table %>%
      gt::tab_header(title = gt::md(table_title))
  }

  gt_table <- gt_table %>%
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
    gt::cols_align(align = "left", columns = "pretty_variable") %>%
    gt::cols_align(align = "center", columns = c("mean_balance", "variance_balance"))

  # Hide appropriate columns
  columns_to_hide <- c("row_id", "poor_balance")
  if (!use_groups) {
    columns_to_hide <- c(columns_to_hide, "shade_flag")
    if (n_models == 1 && is.null(model_labels)) {
      columns_to_hide <- c(columns_to_hide, "model_display")
    }
  }
  if (!has_variance_data) {
    columns_to_hide <- c(columns_to_hide, "variance_balance")
  }

  gt_table <- gt_table %>% gt::cols_hide(columns = columns_to_hide)

  # [Keep all your source note and output code]
  # ...

  return(gt_table)
}
