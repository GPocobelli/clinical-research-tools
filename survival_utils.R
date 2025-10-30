


# ______________________________________________________
#
## Survival Analysis ----
# ______________________________________________________




# -------- required packages ------------

required_setup_packages <- c("survival", "survminer", "tidyverse", "ggpubr", "ggpp")

for (pkg in required_setup_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    install.packages(pkg)
  }
}








#' Calculate time difference between two dates
#'
#' @description
#' Calculate difference in months between two dates
#'
#' @param date1   Date. First date (usually later date).
#' @param date2   Date. Second date (earlier date).
#' @param unit    String. One of "days", "months", or "years". (Default "months").
#'
#' @description   For "years", uses 365.25, for "months", uses 365.25/12
#'
#' @return        Numeric vector with the same length as date1/date2. Difference in months.
#' @export
calc_time <- function(date1, date2, unit = "months") {

  stopifnot(inherits(date1, c("Date", "POSIXct", "POSIXt")),
            inherits(date2, c("Date", "POSIXct", "POSIXt")))

  diff_days <- as.numeric(difftime(date1, date2, units = "days"))
  switch(tolower(unit),
         days   = diff_days,
         months = diff_days / (365.25 / 12),
         years  = diff_days / 365.25,
         stop("`unit` must be 'days', 'months', or 'years'."))
}













#' (Helper) Transform the x-achsis of KM curves
#'
#' @description
#' Multiplication factor to transform the x-axis of KM curves
#' define x-Scale
#'
#' @param xscale   Character of the form "d_m" (= days to months) or "m_y"
#'               (= months to years). See details.
#'
#' @details
#' Allowed pairs (`from`_`to`): d_m, d_y, m_d, m_y, y_d, y_m.
#'
#' @return         Numeric factor
#' @export
get_xscale <- function(xscale = NULL){

  if (is.null(xscale) || xscale %in% c("", "none")) {
    return (1)
  }

  pair <- tolower(xscale)
  if (!grepl("^[dmy]_[dmy]$", pair))
    stop("❗'xscale' must have the pattern: 'd_m', 'm_y', 'm_m', ...")


  # calculations for the right scale
  switch(pair,
         d_m = 12/365.25,
         d_y = 1/365.25,
         m_d = 365.25/12,
         m_y = 1/12,
         y_d = 365.25,
         y_m = 12,
         d_d = 1,
         m_m = 1,
         y_y = 1,
         stop("❗Unknown combination in 'get_xscale()'.")
  )

}












#' (Helper) Define Table theme
#'
#' @param base_size   Integer. Defines the fontsize
#'
#' @returns
#' @export
#'
make_table_theme <- function(base_size = 3) {

  gridExtra::ttheme_minimal(
    core = list(bg_params = list(fill = "white", col = NA),
                fg_params = list(fontface = "plain", fontsize = base_size)),
    colhead = list(bg_params = list(fill = "gray93", col = NA),
                   fg_params = list(fontface = "bold",  fontsize = base_size)),
    rowhead = list(bg_params = list(fill = "gray93", col = NA),
                   fg_params = list(fontface = "plain", fontsize = base_size))
  )
}



















#' (Helper) Median table
#'
#' @description
#' Median Survival (with 95% CI) as a one‑row `tibble`
#'
#' @param fit               `survfit` object.
#' @param xscale            String; Transformation specification passed to `get_xscale()`.
#'                          Character of the form "d_m" (= days to months) or "m_y" (= months to years)
#'                          or "m_m" (= remains months).
#'                          Allowed pairs (`from_to`): d_m, d_y, m_d, m_y, y_d, y_m, d_d, m_m, y_y.
#'                          The description in the table; e.g. "Years", "Months" or "Days".
#'                          Is extracted via xscale automatically.
#'
#'
#' @returns                 A tibble with a single column that is ready to be printed or annotated.
#' @import tidyverse
#' @export
get_median_table <- function(fit, xscale = "m_y") {

  xtrans <- get_xscale(xscale)


  time_unit <- case_when(str_sub(xscale, -1) == "d" ~ "Days",
                         str_sub(xscale, -1) == "m" ~ "Months",
                         str_sub(xscale, -1) == "y" ~ "Years")


  # summary object
  tab <- summary(fit)$table


  # Normalize to data.frame with rows = strata (or 1 row if unstratified)
  if (is.null(dim(tab))) {
    tab <- as.data.frame(t(tab))
    rn  <- "Overall"
  } else {
    tab <- as.data.frame(tab)
    rn  <- rownames(tab)
  }
  # Extract group labels from rownames like "group=Level"
  grp <- sub("^.*?=", "", rn)
  grp[grp == rn] <- rn  # handle "Overall" or already-clean names

  # Compute & format median with CI; handle NR
  med   <- tab[,"median"] * xtrans
  lcl   <- tab[,"0.95LCL"] * xtrans
  ucl   <- tab[,"0.95UCL"] * xtrans

  fmt   <- ifelse(is.finite(med),
                  sprintf("%.2f [%.2f–%.2f]", med, lcl, ucl),
                  "NR")

  out <- tibble::tibble(
    `Group ` = grp,
    Median = fmt
  )

  names(out)[2] <- paste0("Median Survival, \n", time_unit, " [95%-CI]")
  out
}













#' (Helper) Follow-up times table
#'
#' @description
#' Survival probabilities at pre‑defined follow‑up times
#' Follow up Survival probability (in %)
#'
#' @param fit              `survfit` object.
#' @param times            Numeric vector of *original* time scale (e.g. months).
#' @param years            Default; Numeric vector of time scale in years.
#' @param xscale           Transformation specification passed to `get_xscale()`
#'                         The description in the table; e.g. "Years", "Months" or "Days".
#'                         Is extracted via xscale automatically.
#'
#' @returns                A table with two columns (follow up survival probabilities (+ confidence intervals)
#'                         at 1, 2, and 5 years.)
#' @import tibble
#' @export
get_surv_times <- function(fit,
                           times = NULL,
                           years = c(1, 2, 5),
                           xscale = "m_y") {


  xtrans <- get_xscale(xscale)
  origin <- substr(ifelse(is.null(xscale), "y_y", xscale), 1, 1)   # first letter encodes the original unit

  time_unit <- case_when(str_sub(xscale, -1) == "d" ~ "Days",
                         str_sub(xscale, -1) == "m" ~ "Months",
                         str_sub(xscale, -1) == "y" ~ "Years")


  if (is.null(times)) {

    per_year <- switch(origin,
                       d = 365.25,  # days per year
                       m = 12,      # months per year
                       y = 1,       # already in years
                       stop("Cannot determine origin unit from `xscale`."))
    times <- years * per_year
  }


  # summary object
  surv_summary <- summary(fit, times = times, extend = TRUE)


  # Group labels (e.g., "group=A" -> "A")

  grp <- if (!is.null(surv_summary$strata)) {
    paste0(sub("^.*?=", "", surv_summary$strata), ", [95%-CI]")
  } else {
    rep("Overall, [95%-CI]", length(surv_summary$time))
  }


  FU_id  <- surv_summary$time * xtrans
  FU_lab <- format(round(FU_id, 0), trim = TRUE, nsmall = 0)

  surv_txt <- sprintf("%.2f%% [%.2f–%.2f]",
                      surv_summary$surv  * 100,
                      surv_summary$lower * 100,
                      surv_summary$upper * 100)


  df_long <- tibble::tibble(
    `Follow up in {unit}` = FU_lab,
    .FU_id = FU_id,
    Group  = grp,
    Surv   = surv_txt
  )


  # Make a compact wide table: one row per Group, one column per FU point
  df_wide <- tidyr::pivot_wider(
    df_long,
    id_cols    = c(`Follow up in {unit}`, .FU_id),
    names_from = Group,
    values_from= Surv,
    values_fn  = list(Surv = ~ dplyr::first(na.omit(.))),
    values_fill = "_"
  )

  # Ensure stable column order (Group first)
  df_wide <- dplyr::arrange(df_wide, .FU_id)
  df_wide$.FU_id <- NULL
  names(df_wide)[1] <- paste0("Follow up\n in ", time_unit)

  as.data.frame(df_wide, stringsAsFactors = FALSE)
}


















#' helper functions for the cox table

sanitize_id <- function(s) gsub("[^[:alnum:]]+", "", tolower(s))




#' helper
pretty_cox_name <- function(rn, cox_fit) {
  if (grepl("=", rn, fixed = TRUE)) {
    parts <- strsplit(rn, "=", fixed = TRUE)[[1]]
    var <- stringr::str_trim(parts[1]); lvl <- stringr::str_trim(parts[2])
    return(paste0(stringr::str_to_sentence(var), " (", lvl, ")"))
  }
  
  xl <- cox_fit$xlevels
  if (length(xl)) {
    srn <- sanitize_id(rn)
    for (v in names(xl)) {
      for (lvl in xl[[v]]) {
        cand <- sanitize_id(paste0(v, lvl))
        if (identical(cand, srn)) {
          return(paste0(stringr::str_to_sentence(v), " (", lvl, ")"))
        }
      }
    }
  }
  
  if (grepl("[A-Z]", rn)) {
    sp <- gsub("([a-z])([A-Z])", "\\1 \\2", rn)
    parts <- unlist(strsplit(sp, " +"))
    if (length(parts) >= 2) {
      var <- paste(parts[-length(parts)], collapse = " ")
      lvl <- parts[length(parts)]
      return(paste0(stringr::str_to_sentence(var), " (", lvl, ")"))
    }
  }
  stringr::str_to_sentence(rn)
}

















#' @param cox_fit        coxph object
#' @param conf.level     numeric, confidence level for Wald CI (default 0.95)
#' @param digits_hr      int, rounding for HR & CI
#' @param digits_p       int, rounding for p-values (sehr kleine p als "<1e-xx")
#'
#' @return               data.frame with columns: Gruppe, `HR [95% CI]`, p
#' @export
get_cox_table <- function(cox_fit,
                          conf.level = 0.95,
                          digits_hr  = 2,
                          digits_p   = 3) {

  stopifnot(inherits(cox_fit, "coxph"))

  s  <- summary(cox_fit)
  coef <- s$coefficients
  if (is.null(coef) || nrow(coef) == 0L) {
    return(data.frame(Gruppe = character(),
                      `HR [95% CI]` = character(),
                      p = character(),
                      check.names = FALSE))
  }

  # Extract all variables from the cox model output
  zcrit <- stats::qnorm(1 - (1 - conf.level)/2)
  logHR <- coef[, "coef"]
  SE    <- coef[, "se(coef)"]
  lo    <- logHR - zcrit * SE
  hi    <- logHR + zcrit * SE

  HR    <- exp(logHR)
  HR_lo <- exp(lo)
  HR_hi <- exp(hi)

  # Format-helper
  fmt_num <- function(x, d) format(round(x, d), nsmall = d, trim = TRUE, scientific = FALSE)
  fmt_p   <- function(p) {
    if (is.na(p)) return(NA_character_)
    if (p < 10^-(digits_p+1)) {
      paste0("<", format(10^-(digits_p+1), scientific = TRUE))
    } else {
      format(round(p, digits_p), nsmall = digits_p, trim = TRUE, scientific = FALSE)
    }
  }

  term_names <- rownames(coef)
  pretty_term <- function(x) {
    if (grepl("=", x)) {
      sub("^([^=]+)=(.*)$", "\\1: \\2", x)
    } else {
      x
    }
  }
  grp <- vapply(term_names, pretty_cox_name, character(1), cox_fit = cox_fit)
  
  grp <- stringr::str_to_sentence(grp)

  hr_txt <- paste0(fmt_num(HR, digits_hr), " [",
                   fmt_num(HR_lo, digits_hr), "–",
                   fmt_num(HR_hi, digits_hr), "]")

  p_txt  <- vapply(coef[, "Pr(>|z|)"], fmt_p, character(1))

  df <- data.frame(
    Strata = grp,
    check.names = FALSE,
    `HR [95% CI]` = hr_txt,
    p = p_txt,
    stringsAsFactors = FALSE
  )

  rownames(df) <- NULL
  df
}












#' (Helper) Normalize show_tbls specification
#'
#' @param show_tbls     various types: logical scalar; logical vector (len 2/3);
#'                      named logical vector with names among c("median","followup","cox");
#'                      character vector with any of those names.
#' @param has_cox       logical; whether a Cox fit was provided
#' @return              named logical vector c(median=., followup=., cox=.)
#' @export
normalize_show_tbls <- function(show_tbls, has_cox) {

  valid_names <- c("median", "followup", "cox")

  # defaults mimic previous behavior: median & followup shown, cox shown if available
  defaults <- c(median = TRUE, followup = TRUE, cox = has_cox)

  #  scalar logical
  if (is.logical(show_tbls) && length(show_tbls) == 1L) {
    return(setNames(rep(isTRUE(show_tbls), 3), valid_names))
  }

  # character vector (names of tables to show)
  if (is.character(show_tbls)) {
    nm <- tolower(show_tbls)
    if (!all(nm %in% valid_names)) {
      stop("`show_tbls` character values must be among: ", paste(valid_names, collapse = ", "))
    }
    out <- setNames(rep(FALSE, 3), valid_names)
    out[nm] <- TRUE
    return(out)
  }

  # logical vector (possibly named) length 2 or 3
  if (is.logical(show_tbls)) {
    # named logical: use provided names; others keep default
    if (!is.null(names(show_tbls))) {
      nm <- tolower(names(show_tbls))
      if (!all(nm %in% valid_names)) {
        stop("`show_tbls` named logical must use names among: ", paste(valid_names, collapse = ", "))
      }
      out <- defaults
      out[nm] <- as.logical(show_tbls)
      return(out)
    }

    # unnamed: length 2 -> c(median, followup), cox = default (has_cox)
    if (length(show_tbls) == 2L) {
      return(c(median  = as.logical(show_tbls[1]),
               followup= as.logical(show_tbls[2]),
               cox     = has_cox))
    }

    # unnamed: length 3 -> c(median, followup, cox)
    if (length(show_tbls) == 3L) {
      return(setNames(as.logical(show_tbls), valid_names))
    }

    stop("`show_tbls` logical must be length 1, 2, or 3 (optionally named).")
  }

  stop("Unsupported `show_tbls` type. Use TRUE/FALSE, a logical vector (len 2/3), a named logical, or a character vector of names.")
}

















#' High‑level Kaplan–Meier plot wrapper
#'
#' @description
#' High‑level Kaplan–Meier plot wrapper with add-ons like median & follow-up times table,
#' table positions, customazation of the scale, title, title-size and more.
#'
#' @param data                     Optional `data.frame` containing all variables. used for the risk table. If omitted the
#'                                 risk table is suppressed (because `ggsurvplot()` requires `data` for that panel).
#' @param fit                      `survfit` object. Stratification present in `fit` will be respected automatically.
#' @param time_unit                String. Original time unit of the model (informational, used for axis label).
#' @param xscale                   See `get_xscale()`. Default = "m_y" (= from months to years).
#'                                 Allowed pairs (`from_to`): "d_m", "d_y", "m_d", "m_y", "y_d", "y_m", "d_d", "m_m", "y_y",
#'                                 or "" / NULL (= no conversion).
#' @param x_end                    Numeric. Optional. Gets calculated by default from `survfit` object.
#'                                 Defines the x-scale end point. In *original* time unit.
#' @param scale_break              Numeric. Break points of the x-scale. Default = 24 months. In *original* time unit.
#' @param title                    String. Plot title
#' @param title_size               Numeric. Font size of the title. Default = 11.
#' @param show_tbls                Controls which summary tables are drawn inside the KM-plot.
#'                                  Accepts:
#'                                   - `TRUE` / `FALSE`: show all / none (including Cox if `cox_fit` is provided)
#'                                   - unnamed logical length 2: `c(median, followup)` (Cox uses default = shown if available)
#'                                   - unnamed logical length 3: `c(median, followup, cox)`
#'                                   - named logical: any of `c(median = TRUE, followup = FALSE, cox = TRUE)`
#'                                   - character vector: any of `c("median", "followup", "cox")`
#'                                  Examples:
#'                                     `show_tbls = TRUE`; `show_tbls = c(FALSE, TRUE)`;
#'                                     `show_tbls = c(median = TRUE, followup = FALSE, cox = TRUE)`;
#'                                     `show_tbls = c("median","cox")`.
#' @param tbl1_pos                 Numeric Vector length 2 (x, y) for `get_median_table()`.
#' @param tbl2_pos                 Numeric Vector length 2 (x, y) for `get_surv_times()`.
#' @param followup_times           Numeric vector of individual length for `get_surv_times()`. Defaut = NULL.
#' @param risk_table_size          Numeric. Font size of the risk table.
#' @param palette                  Vector. Arguments for color-pallett group-specific.
#' @param cox_fit                  `coxph` object.
#' @param cox_tbl_pos              Numeric Vector length 2 (x, y). Default = c(0.9, 0.75)
#' @param ...                      Further arguments forwarded to `ggsurvplot()`
#'                                 (e.g. `conf.int = TRUE`, `pval = TRUE`, `linetype = "strata"`).
#'
#' @returns                        A `ggsurvplot` object list.
#' @export
#' @import survminer
#' @import tidyverse
#' @import survival
#' @import ggpp
#' @import ggpubr
#' @examples
#' library(survival)
#' surv <- Surv(time = lung$time, event = lung$status)
#' fit <- survfit(surv ~ sex, data = lung)
#' cox <- coxph(surv ~ sex, data = lung)
#' create_surv_plot(
#'    data = lung,
#'    fit = fit,
#'    xscale = "d_m",              # days to months
#'    time_unit = "Months",
#'    scale_break = 182.625,       # 6 months (182.625 days)
#'    x_end = 1000,                # optional: end of the x-axis in days
#'    title = "Kaplan–Meier Curve",
#'    tbl2_pos = c(0.9, 1.2),     # position of the follow up table
#'    risk_table_size = 4,         # fontsize of the values
#'    cox_fit = cox,
#'    cox_tbl_pos = c(0.5, 0),
#'    show_cox_if_grouped = T,
#'    show_tbls = c(T, TRUE)
#' )

create_surv_plot <- function(data        = NULL,
                             fit,
                             xscale      = "m_y",
                             scale_break = 24,
                             x_end       = NULL,
                             title       = "Kaplan-Meier curve",
                             title_size  = 11,
                             show_tbls   = TRUE,
                             tbl1_pos    = c(0, 0),
                             tbl2_pos    = c(0.8, 0),
                             followup_times  = NULL,
                             risk_table_size = 3.5,
                             palette     = NULL,
                             cox_fit     = NULL,
                             cox_tbl_pos = c(0.9, 0.75),
                             ...) {

  stopifnot(inherits(fit, "survfit"))


  x_end_raw     <- if (is.null(x_end)) dplyr::last(fit$time) else x_end
  x_end_aligned <- scale_break * ceiling(x_end_raw / scale_break)


  reldata_x <- function(rel) rel * x_end_aligned
  reldata_y <- function(rel) rel



  time_unit <- case_when(str_sub(xscale, -1) == "d" ~ "Days",
                         str_sub(xscale, -1) == "m" ~ "Months",
                         str_sub(xscale, -1) == "y" ~ "Years")



  risk_tbl <- if (!is.null(data)) "nrisk_cumevents" else FALSE

  dots <- list(...)

  base_args <- list(
    fit        = fit,
    data       = data,
    risk.table = risk_tbl,
    xlab       = paste0("Time (", time_unit, ")"),
    ylab       = "Survival probability (in %)",
    title      = title,
    break.time.by       = scale_break,
    risk.table.fontsize = risk_table_size,
    censor.shape        = "|",
    ggtheme             = ggpubr::theme_pubr(),
    palette             = palette
  )


  has_strata <- !is.null(fit$strata) && length(fit$strata) > 0

  n_groups <- if (is.null(fit$strata)) 1L else length(fit$strata)

  if (n_groups == 1L) {
    # Only set if the user has not specified otherwise
    if (is.null(dots$palette)) base_args$palette <- "black"
    if (is.null(dots$pval))    base_args$pval    <- FALSE
  }


  if (has_strata) {
    labs_raw   <- names(fit$strata)
    labs_clean <- gsub("(^|, )[^=]+=", "", labs_raw)
    if (is.null(dots$legend.labs)) base_args$legend.labs   <- labs_clean
    if (is.null(dots$legend.title)) base_args$legend.title <- ""
  }

  if (is.null(dots$surv.median.line))       base_args$surv.median.line       <- "hv"
  if (is.null(dots$conf.int))               base_args$conf.int               <- TRUE
  if (has_strata && is.null(dots$pval))     base_args$pval                   <- TRUE
  if (is.null(dots$risk.table))             base_args$risk.table             <- risk_tbl
  if (is.null(dots$risk.table.col))         base_args$risk.table.col         <- if (has_strata) "strata" else "black"

  if (is.null(dots$risk.table.y.text))      base_args$risk.table.y.text      <- TRUE
                                       else base_args$risk.table.y.text      <- isTRUE(dots$risk.table.y.text)

  if (is.null(dots$risk.table.y.text.col))  base_args$risk.table.y.text.col  <- TRUE
                                       else base_args$risk.table.y.text.col  <- isTRUE(dots$risk.table.y.text.col)

  if (is.null(dots$xlim))                   base_args$xlim                   <- c(0, x_end_aligned)


  # Merge user args with `ggsurvplot` args
  args <- utils::modifyList(base_args, dots, keep.null = TRUE)


  # Axis transformation
  breaks_x <- seq(0, x_end_aligned, by = scale_break)
  xtrans   <- get_xscale(xscale)
  breaks_y <- seq(0, 1.5, by = 0.25)

  plot_obj <- do.call(survminer::ggsurvplot, args)


  y_headroom <- 0.25

  # Scale Adjustments
  plot_obj$plot <- plot_obj$plot +
    ggplot2::scale_x_continuous(breaks = breaks_x,
                                labels = round(breaks_x * xtrans, 2),
                                expand = c(0.05, 0.2)) +

    ggplot2::scale_y_continuous(limits = c(0, 1 + y_headroom - 0.01),
                                labels = function(x) x * 100,
                                expand = ggplot2::expansion(mult = c(0.01, 0))) +

    ggplot2::coord_cartesian(xlim = c(0, x_end_aligned), clip = "on") +

    ggplot2::theme(axis.text.x  = element_text(size = title_size),
                   axis.title.x = element_text(size = title_size, face = "bold"),
                   axis.text.y  = element_text(size = title_size),
                   axis.title.y = element_text(size = title_size, face = "bold"),
                   title        = element_text(size = title_size, face = "bold"),
                   legend.position = if (has_strata && n_groups > 1L) "bottom" else "none",
                   legend.direction	= "horizontal",
                   legend.justification = c(-0.01, 0.5),
                   legend.box       = "horizontal"
    ) +
    ggplot2::guides(color = guide_legend(nrow = 2, byrow = TRUE))







  # ---- create tables ----
  has_cox <- !is.null(cox_fit)
  show_vec <- normalize_show_tbls(show_tbls, has_cox = has_cox)
  show_median   <- isTRUE(show_vec["median"])
  show_followup <- isTRUE(show_vec["followup"])
  show_cox      <- isTRUE(show_vec["cox"])


  ## ---- Create & annotate Median / Follow-up tables ----
  if (show_median || show_followup) {

    if (show_median) {
      tbl1 <- as.data.frame(get_median_table(fit, xscale = xscale))
      plot_obj$plot <- plot_obj$plot +
        ggpp::annotate(geom = "table",
                       x = reldata_x(tbl1_pos[1]),
                       y = reldata_y(tbl1_pos[2]),
                       label = tbl1,
                       table.theme = make_table_theme(title_size - 2))
    }

    if (show_followup) {
      tbl2 <- as.data.frame(get_surv_times(fit, times = followup_times, xscale = xscale))
      plot_obj$plot <- plot_obj$plot +
        ggpp::annotate(geom = "table",
                       x = reldata_x(tbl2_pos[1]),
                       y = reldata_y(tbl2_pos[2]),
                       label = tbl2,
                       table.theme = make_table_theme(title_size - 2))
    }
  }


  ## ---- Cox-model table ----
  if (show_cox && has_cox && has_strata) {
    ctab <- as.data.frame(get_cox_table(cox_fit))
    plot_obj$plot <- plot_obj$plot +
      ggpp::annotate(
        geom = "table",
        x = reldata_x(cox_tbl_pos[1]), y = reldata_y(cox_tbl_pos[2]),
        label = ctab,
        table.theme = make_table_theme(title_size - 2)
      )
  }





# ---- risk table scales ----
  if (!isFALSE(args$risk.table) && !is.null(plot_obj$table)) {


    plot_obj$table <- plot_obj$table +
      ggplot2::scale_x_continuous(
        limits = c(0, x_end_aligned),
        breaks = breaks_x,
        labels = round(breaks_x * xtrans, 2),
        expand = c(0.05, 0)
      ) +
      ggplot2::xlab("") +
      ggplot2::theme(
        axis.title.y = ggplot2::element_blank(),
        axis.ticks.y = ggplot2::element_blank(),
        axis.ticks.x = ggplot2::element_blank(),
        legend.position = "none"
      )

    if (!isFALSE(args$risk.table) &&
        isFALSE(args$risk.table.y.text) &&
        isTRUE(args$risk.table.y.text.col) &&
        !is.null(plot_obj$table)) {

      ## ---- Strata levels ----
      all_layers <- lapply(seq_along(plot_obj$table$layers),
                           function(i) ggplot2::layer_data(plot_obj$table, i))
      strata_vec <- unlist(lapply(all_layers, function(d)
        if ("strata" %in% names(d)) as.character(d$strata) else NULL),
        use.names = FALSE)

      strata_levels <- unique(stats::na.omit(strata_vec))
      strata_levels <- strata_levels[nchar(strata_levels) > 0]

      # just apply if there are more than 2 groups
      if (length(strata_levels) >= 2) {

        # order like before
        strata_levels <- rev(strata_levels)

        # Adjust space left of the table, without changeing something of the y-achsis
        left_pad <- 0.04 * x_end_aligned   # 4% of the width
        plot_obj$table <- plot_obj$table +
          ggplot2::coord_cartesian(xlim = c(-left_pad, x_end_aligned), clip = "off")


        marker_df <- data.frame(
          .y = factor(strata_levels, levels = strata_levels),
          .x = -left_pad * 0.5
        )

        # ---- colored label of each strata ----
        plot_obj$table <- plot_obj$table +
          ggplot2::geom_tile(
            data = marker_df,
            mapping = ggplot2::aes(x = .x, y = .y, fill = .y),
            width = left_pad * 0.8,
            height = 0.6,
            color = NA,
            inherit.aes = FALSE,
            show.legend = FALSE
          )

        # ---- same palette from the plot ----
        # same colors as the curves
        if (!is.null(args$palette)) {
          if (is.atomic(args$palette)) {
            vals <- rep_len(args$palette, length(strata_levels))
            plot_obj$table <- plot_obj$table +
              ggplot2::scale_fill_manual(values = vals, guide = "none")
          }
        }
      }
    }
  }
  return(plot_obj)
}


























