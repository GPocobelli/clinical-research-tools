


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
# get_xscale <- function(xscale = "m_y"){
#   # calculations for the right scale
#   xtrans <- switch(xscale,
#                    d_m = 12/365.25,
#                    d_y = 1/365.25,
#                   m_d = 365.25/12,
#                    m_y = 1/12,
#                    y_d = 365.25,
#                    y_m = 12,
#                    1
#   )
#   return(xtrans)
# }



# (Helper) Transform the x-achsis of KM curves
#
# @description
# Multiplication factor to transform the x-axis of KM curves
# define x-Scale
#
# @param xscale   Character of the form "d_m" (= days to months) or "m_y"
#                 (= months to years) or "m_m" (= remains months). See details.
#
# @details
# Allowed pairs (`from_to`): "d_m", "d_y", "m_d", "m_y", "y_d", "y_m", "d_d", "m_m", "y_y",
# or "" / NULL (= no conversion).
#
# @return         Numeric factor (1 = no conversion)
# @export
# get_xscale <- function(xscale = "m_y"){
#   # calculations for the right scale
#   xtrans <- switch(xscale,
#                    d_m = 12/365.25,
#                    d_y = 1/365.25,
#                    m_d = 365.25/12,
#                    m_y = 1/12,
#                    y_d = 365.25,
#                    y_m = 12,
#                    1
#   )
#   return(xtrans)
# }











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
    Group = grp,
    Median = fmt
  )

  names(out)[2] <- paste0("Median Survival, ", time_unit, "\n[95%-CI]")
  out


  # # Get values
  # records <- tab["records"]
  # events  <- tab["events"]
  # median  <- round(tab["median"] * xtrans, 2)
  # lcl     <- round(tab["0.95LCL"] * xtrans, 2)
  # ucl     <- round(tab["0.95UCL"] * xtrans, 2)
  #
  # # Output
  # data.frame(
  #   Median_survival = paste0(median, "  [", lcl, "–", ucl, "]"),
  #   check.names = FALSE  # prevents conversion of column name
  # ) %>%
  #   setNames(paste0("Median\n Survival, ", time_unit,"\n[95%-CI]"))
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
    sub("^.*?=", "", surv_summary$strata)
    }
  else {
    rep("Overall", length(surv_summary$time))
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
  names(df_wide)[1] <- paste0("Follow up in ", time_unit)

  as.data.frame(df_wide, stringsAsFactors = FALSE)



  # ---

  # summary object surv_summary <- summary(fit, times = times) # Group labels (e.g., "group=A" -> "A") if (!is.null(surv_summary$strata)) { grp <- sub("^.*?=", "", surv_summary$strata) } else { grp <- rep("Overall", length(surv_summary$time)) } df_long <- tibble::tibble( Group = grp, FU = round(surv_summary$time * xtrans, 2), Surv = sprintf("%.2f%% [%.2f–%.2f]", surv_summary$surv * 100, surv_summary$lower * 100, surv_summary$upper * 100) ) # Make a compact wide table: one row per Group, one column per FU point df_wide <- tidyr::pivot_wider( df_long, id_cols = Group, names_from = FU, values_from = Surv, names_sort = TRUE, names_glue = paste0("Follow up = ", "{FU}", " ", "\nSurvival (%), [95%-CI]") ) # Ensure stable column order (Group first) df_wide <- df_wide[, c("Group", setdiff(names(df_wide), "Group"))] df_wide



  #
  # time_months <- round(surv_summary$time, 2)
  # surv_percent <- round(surv_summary$surv * 100, 2)
  # ci_lower <- round(surv_summary$lower * 100, 2)
  # ci_upper <- round(surv_summary$upper * 100, 2)
  #
  #
  # # table
  # tibble::tibble(
  #   FU_in_Years = round(time_months * xtrans, 2),
  #   Survival =  paste0(surv_percent, "% [", ci_lower, "–", ci_upper, "]")) %>%
  #
  #   rename_with(~ c(paste0("Follow up\nin ", time_unit),
  #                   paste0("Survival probability (%),\n[95%-CI]")))
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
#' @param title_size               Numeric. Font size of the title. Defaule = 2.7.
#' @param show_tbls                Logical. Adds summary tables inside the Kaplan-Meier-plot.
#' @param tbl1_pos                 Numeric Vector length 2 (x, y) for `get_median_table()`.
#' @param tbl2_pos                 Numeric Vector length 2 (x, y) for `get_surv_times()`.
#' @param followup_times           Numeric vector of individual length for `get_surv_times()`. Defaut = NULL.
#' @param risk_table_size          Numeric. Font size of the risk table.
#' @param palette                  Vector. Arguments for color-pallett group-specific.
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
#' fit <- survfit(surv ~ 1, data = lung)
#' create_surv_plot(
#'    data = lung,
#'    fit = fit,
#'    xscale = "d_m",              # days to months
#'    time_unit = "Months",
#'    scale_break = 182.625,       # 6 months (182.625 days)
#'    x_end = 1000,                # optional: end of the x-axis in days
#'    title = "Kaplan–Meier Curve",
#'    tbl2_pos = c(0.8, 0.99),     # position of the follow up table
#'    risk_table_size = 4          # fontsize of the values
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
                             show_cox_if_grouped = TRUE,
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
    # nur setzen, wenn der/die Nutzer:in nichts anderes vorgeben hat
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


  # Nutzer-Argumente überschreiben Defaults (und vermeiden doppelte Namen)
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

    ggplot2::coord_cartesian(xlim = c(0, x_end_aligned), clip = "off") +

    ggplot2::theme(axis.text.x  = element_text(size = title_size),
                   axis.title.x = element_text(size = title_size, face = "bold"),
                   axis.text.y  = element_text(size = title_size, face = "bold"),
                   axis.title.y = element_text(size = title_size, face = "bold"),
                   title        = element_text(size = title_size, face = "bold"),
                   legend.position = if (has_strata && n_groups > 1L) "right" else "none",
                   legend.justification = c(1, 0)
    ) +
    ggplot2::guides(color = guide_legend(nrow = 3, ncol = 2))


  # create tables
  if (show_tbls){
    # get the Median time
    tbl1 <- get_median_table(fit, xscale = xscale)
    tbl1 <- as.data.frame(tbl1)
    # get the Followup Survival Probability in %
    tbl2 <- get_surv_times(fit, times = followup_times, xscale = xscale)
    tbl2 <- as.data.frame(tbl2)


    # annotate them in the plot
    plot_obj$plot <- plot_obj$plot +
      ggpp::annotate(geom = "table",
                     x = reldata_x(tbl1_pos[1]),
                     y = reldata_y(tbl1_pos[2]),
                     label = tbl1,
                     table.theme = make_table_theme(title_size - 2)) +

      ggpp::annotate(geom = "table",
                     x = reldata_x(tbl2_pos[1]),
                     y = reldata_y(tbl2_pos[2]),
                     label = tbl2,
                     table.theme = make_table_theme(title_size - 2))
  }


  # Cox-model table
  if (!is.null(cox_fit) && (has_strata || !show_cox_if_grouped)) {
    ctab <- as.data.frame(get_cox_table(cox_fit))
    plot_obj$plot <- plot_obj$plot +
      ggpp::annotate(
        geom = "table",
        x = reldata_x(cox_tbl_pos[1]), y = reldata_y(cox_tbl_pos[2]),
        label = ctab, table.theme = make_table_theme(title_size - 2))
  }


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

      all_layers <- lapply(seq_along(plot_obj$table$layers),
                           function(i) ggplot2::layer_data(plot_obj$table, i))
      strata_vec <- unlist(lapply(all_layers, function(d)
        if ("strata" %in% names(d)) as.character(d$strata) else NULL),
        use.names = FALSE)

      strata_levels <- unique(stats::na.omit(strata_vec))
      strata_levels <- strata_levels[nchar(strata_levels) > 0]

      if (length(strata_levels) >= 2) {
        strata_levels <- rev(strata_levels)
        left_pad <- 0.04 * x_end_aligned   
        plot_obj$table <- plot_obj$table +
          ggplot2::coord_cartesian(xlim = c(-left_pad, x_end_aligned), clip = "off")

        marker_df <- data.frame(
          .y = factor(strata_levels, levels = strata_levels),
          .x = -left_pad * 0.5
        )

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













