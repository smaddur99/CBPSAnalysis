#' Export NPCBPS / CBPS Results to an Excel Summary Table
#'
#' Collects the pooled results and covariate balance from one or more
#' \code{npcbps_weighted_analysis()} or \code{cbps_weighted_analysis()} runs and writes
#' them to an .xlsx file with the same layout as the hand-built results table:
#' one row per model term, with model information and balance on the first row
#' of each model.
#'
#' Built on \code{extract_cbps_results()}: each results object is passed through it,
#' so the Excel table and the R summaries come from the same numbers. You can also
#' pass objects that were already run through \code{extract_cbps_results()}.
#'
#' Columns (in order):
#' \itemize{
#'   \item Hypothesis, model name, treatment, outcome, sens
#'   \item variable, estimate, se, ci low, ci high, p-value, n: pooled (Rubin's rules)
#'     results for each non-intercept term, treatment first
#'   \item For each weighting covariate k: covar. k, smd/ R k, vr min k, vr max k,
#'     vr balanced k. smd/ R is |SMD| (binary treatment) or |r| (continuous treatment)
#'     after weighting, across imputations (worst case by default). Variance ratios are
#'     the min and max across imputations; "N/A" for continuous treatments and binary covariates.
#'   \item Optional extras at the end: family, coef scale, exp estimate, exp ci low,
#'     exp ci high (odds ratios for ordinal / logit models), se method, balance stat
#' }
#'
#' @param results A results object, or a list of them. List names become the
#'   "model name" column (e.g. \code{list(preg_advice = res_adv, mpss = res_mpss)}).
#'   Each element can be raw output from \code{npcbps_weighted_analysis()} /
#'   \code{cbps_weighted_analysis()} or output from \code{extract_cbps_results()}.
#' @param file Character. Path of the .xlsx file to write.
#' @param hypothesis Character vector, one per model (or length 1, recycled). Default NA.
#' @param sens Character vector, one per model (or length 1, recycled). Default "main".
#' @param treatment,outcome Optional character vectors, one per model. By default they
#'   are read from the fitted model's formula (needed for already-extracted inputs).
#' @param covariate_labels Optional named character vector to relabel covariates,
#'   e.g. \code{c(gma_age_t1 = "MGM age", m_age_t1 = "M age")}.
#' @param balance_stat "max" (default; worst imputation) or "mean" across imputations.
#' @param vr_range Numeric length 2. Acceptable variance ratio range (default c(0.5, 2)).
#' @param digits Integer. Significant digits for estimates and variance ratios;
#'   decimal places for balance statistics (default 3).
#' @param small_threshold Numeric. Balance statistics below this are written as
#'   "< 0.001" (default 0.001).
#' @param include_extra Logical. Append the family / odds ratio / SE method columns (default TRUE).
#' @param sheet_name Character. Worksheet name (default "Results").
#' @param overwrite Logical. Overwrite \code{file} if it exists (default TRUE).
#' @param ... Passed to \code{extract_cbps_results()} (e.g. \code{smd_threshold}).
#'
#' @return Invisibly, the table as a data frame (all columns character).
#'
#' @examples
#' \dontrun{
#' export_cbps_excel(
#'   results = list(preg_advice = res_adv, mpss = res_mpss, food_pro_ex_dist = res_food),
#'   file = "npcbps_results.xlsx",
#'   hypothesis = "h1",
#'   sens = "main",
#'   covariate_labels = c(gma_age_t1 = "MGM age", m_age_t1 = "M age", ga_t1_yr_c = "g age")
#' )
#' }
#'
#' @export
export_cbps_excel <- function(results,
                              file,
                              hypothesis = NA_character_,
                              sens = "main",
                              treatment = NULL,
                              outcome = NULL,
                              covariate_labels = NULL,
                              balance_stat = c("max", "mean"),
                              vr_range = c(0.5, 2),
                              digits = 3,
                              small_threshold = 0.001,
                              include_extra = TRUE,
                              sheet_name = "Results",
                              overwrite = TRUE,
                              ...) {

  if (!requireNamespace("openxlsx", quietly = TRUE)) {
    stop("Package openxlsx is required: install.packages(\"openxlsx\")")
  }
  if (!exists("extract_cbps_results", mode = "function")) {
    stop("extract_cbps_results() must be loaded first.")
  }
  balance_stat <- match.arg(balance_stat)

  # ---- one results object or a list of them ----
  is_single <- !is.null(results$pooled_results) || !is.null(results$pooled_estimates)
  if (is_single) results <- list(results)
  n_models <- length(results)

  recycle <- function(x, what) {
    if (is.null(x)) return(rep(NA_character_, n_models))
    if (length(x) == 1) return(rep(x, n_models))
    if (length(x) != n_models) stop("`", what, "` must have length 1 or ", n_models, ".")
    x
  }
  hypothesis <- recycle(hypothesis, "hypothesis")
  sens       <- recycle(sens, "sens")
  treatment  <- recycle(treatment, "treatment")
  outcome    <- recycle(outcome, "outcome")

  model_names <- names(results)
  if (is.null(model_names)) model_names <- rep("", n_models)

  # ---- formatting helpers ----
  is_missing <- function(x) is.null(x) || length(x) == 0 || is.na(x[1])
  fmt_num <- function(x) if (is_missing(x)) "N/A" else signif(x[1], digits)
  fmt_bal <- function(x) {
    if (is_missing(x)) return("N/A")
    if (abs(x[1]) < small_threshold) return(paste("<", small_threshold))
    round(x[1], digits)
  }
  relabel <- function(v) {
    if (!is.null(covariate_labels) && v %in% names(covariate_labels)) covariate_labels[[v]] else v
  }

  rows <- list()
  max_covariates <- 0

  for (i in seq_len(n_models)) {
    res <- results[[i]]

    # ---- run through extract_cbps_results() unless already extracted ----
    is_extracted <- !is.null(res$pooled_estimates) && is.null(res$imputation_results)
    ext <- if (is_extracted) {
      res
    } else {
      extract_cbps_results(res, include_balance_details = TRUE, vr_range = vr_range, ...)
    }

    pe <- ext$pooled_estimates
    if (is.null(pe) || nrow(pe) == 0) {
      warning("Model ", i, " has no pooled estimates; skipped.", call. = FALSE)
      next
    }

    # ---- treatment and outcome names from the fitted model's formula ----
    if (!is_extracted && !is.null(res$model)) {
      tt <- stats::terms(res$model)
      if (is.na(outcome[i]))   outcome[i]   <- all.vars(tt)[1]
      if (is.na(treatment[i])) treatment[i] <- attr(tt, "term.labels")[1]
    }

    pe <- pe[pe$variable != "(Intercept)", , drop = FALSE]
    if (is.na(treatment[i])) treatment[i] <- pe$variable[1]
    is_treat <- startsWith(pe$variable, treatment[i])
    pe <- rbind(pe[is_treat, , drop = FALSE], pe[!is_treat, , drop = FALSE])  # treatment first

    if (model_names[i] == "") model_names[i] <- outcome[i]

    n_final <- ext$sample_sizes$value[ext$sample_sizes$metric == "final_n"]

    # ---- covariate balance, in propensity-formula order ----
    bal_cells <- list()
    bbi <- ext$balance_by_imputation
    has_bal <- !is.null(bbi) && nrow(bbi) > 0 && "variable" %in% names(ext$balance) &&
      "max_balance_stat" %in% names(ext$balance)

    if (has_bal) {
      cov_order <- unique(bbi$variable[bbi$imputation == min(bbi$imputation)])
      bal <- ext$balance[match(cov_order, ext$balance$variable), , drop = FALSE]
      bal <- bal[!is.na(bal$variable) & !(bal$type %in% "Distance"), , drop = FALSE]
      is_binary_treat <- identical(ext$treatment_type, "binary")

      for (k in seq_len(nrow(bal))) {
        stat_k <- if (balance_stat == "max") bal$max_balance_stat[k] else bal$mean_balance_stat[k]
        vr_min <- if (is_binary_treat) bal$min_var_ratio[k] else NA_real_
        vr_max <- if (is_binary_treat) bal$max_var_ratio[k] else NA_real_
        bal_cells[[paste("covar.", k)]]      <- relabel(bal$variable[k])
        bal_cells[[paste("smd/ R", k)]]      <- fmt_bal(stat_k)
        bal_cells[[paste("vr min", k)]]      <- fmt_num(vr_min)
        bal_cells[[paste("vr max", k)]]      <- fmt_num(vr_max)
        bal_cells[[paste("vr balanced", k)]] <- if (is_missing(vr_min)) "N/A" else
          (vr_min >= vr_range[1] & vr_max <= vr_range[2])
      }
      max_covariates <- max(max_covariates, nrow(bal))
    }

    # ---- one row per term; model info and balance on the first row ----
    for (j in seq_len(nrow(pe))) {
      row <- list(
        "variable" = pe$variable[j],
        "estimate" = fmt_num(pe$pooled_estimate[j]),
        "se"       = fmt_num(pe$pooled_se[j]),
        "ci low"   = fmt_num(pe$pooled_conf_low[j]),
        "ci high"  = fmt_num(pe$pooled_conf_high[j]),
        "p-value"  = fmt_num(pe$pooled_p_value[j])
      )
      if (include_extra) {
        row[["exp estimate"]] <- if ("exp_estimate" %in% names(pe)) fmt_num(pe$exp_estimate[j]) else NULL
        row[["exp ci low"]]   <- if ("exp_ci_lower" %in% names(pe)) fmt_num(pe$exp_ci_lower[j]) else NULL
        row[["exp ci high"]]  <- if ("exp_ci_upper" %in% names(pe)) fmt_num(pe$exp_ci_upper[j]) else NULL
        row[["se method"]]    <- if ("se_method" %in% names(pe)) pe$se_method[j] else NULL
      }
      if (j == 1) {
        row <- c(list(
          "Hypothesis" = hypothesis[i],
          "model name" = model_names[i],
          "treatment"  = treatment[i],
          "outcome"    = outcome[i],
          "sens"       = sens[i],
          "n"          = n_final
        ), row, bal_cells)
        if (include_extra) {
          row[["family"]]       <- ext$family_label
          row[["coef scale"]]   <- ext$coef_scale
          row[["balance stat"]] <- ext$balance_statistic
        }
      }
      rows[[length(rows) + 1]] <- row
    }
  }

  if (length(rows) == 0) stop("No results to write.")

  # ---- header ----
  header <- c("Hypothesis", "model name", "treatment", "outcome", "sens",
              "variable", "estimate", "se", "ci low", "ci high", "p-value", "n")
  for (k in seq_len(max_covariates)) {
    header <- c(header, paste("covar.", k), paste("smd/ R", k), paste("vr min", k),
                paste("vr max", k), paste("vr balanced", k))
  }
  if (include_extra) {
    header <- c(header, "family", "coef scale", "exp estimate", "exp ci low",
                "exp ci high", "se method", "balance stat")
  }

  # Covariate slots a model doesn't use are written as "N/A" on its first row
  for (r in seq_along(rows)) {
    if (!is.null(rows[[r]][["Hypothesis"]]) || !is.null(rows[[r]][["model name"]])) {
      for (k in seq_len(max_covariates)) {
        for (h in c(paste("covar.", k), paste("smd/ R", k), paste("vr min", k),
                    paste("vr max", k), paste("vr balanced", k))) {
          if (is.null(rows[[r]][[h]])) rows[[r]][[h]] <- "N/A"
        }
      }
    }
  }

  # ---- write: one cell at a time so numbers stay numbers and "N/A" / "< 0.001" stay text ----
  wb <- openxlsx::createWorkbook()
  openxlsx::modifyBaseFont(wb, fontSize = 12, fontName = "Aptos Narrow")
  openxlsx::addWorksheet(wb, sheet_name)
  openxlsx::writeData(wb, sheet_name, x = matrix(header, nrow = 1),
                      startRow = 1, startCol = 1, colNames = FALSE)

  for (r in seq_along(rows)) {
    for (h in names(rows[[r]])) {
      v <- rows[[r]][[h]]
      col <- match(h, header)
      if (is_missing(v) || is.na(col)) next
      openxlsx::writeData(wb, sheet_name, x = v, startRow = r + 1, startCol = col, colNames = FALSE)
    }
  }

  openxlsx::freezePane(wb, sheet_name, firstRow = TRUE)
  openxlsx::setColWidths(wb, sheet_name, cols = seq_along(header), widths = "auto")
  openxlsx::saveWorkbook(wb, file, overwrite = overwrite)

  # ---- return the same table as a data frame ----
  table_df <- as.data.frame(
    do.call(rbind, lapply(rows, function(row) {
      vapply(header, function(h) {
        v <- row[[h]]
        if (is_missing(v)) NA_character_ else as.character(v)
      }, character(1))
    })),
    stringsAsFactors = FALSE
  )
  names(table_df) <- header

  message("Wrote ", nrow(table_df), " rows (", n_models, " models) to ", file)
  invisible(table_df)
}
