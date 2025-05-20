

### Survival Analysis ----




#' Calculate difference in months between two dates
#'
#' @param date1   Date. First date (usually later date).
#' @param date2   Date. Second date (earlier date).
#'
#' @return        Numeric. Difference in months.
#'
calc_time <- function(date1, date2) {
  as.numeric(round(difftime(date1, date2, units = "days") / (365.25 / 12), digits = 2))
}













#' Get the multiplication factor to convert x-axis for survival plots
#' define x-Scale
#' 
#' @param xscale   Character. Abbreviation for axis transformation, e.g. "m_y" (month to year).
#'
#' @return         Numeric. Multiplication factor for axis.
#'
get_xscale <- function(xscale = "m_y"){
  # calculations for the right scale
  xtrans <- switch(xscale,
                   d_m = 12/365.25,
                   d_y = 1/365.25,
                   m_d = 365.25/12,
                   m_y = 1/12,
                   y_d = 365.25,
                   y_m = 12,
                   1
  )
  return(xtrans)
}














# Get the Median Survival Table
#' Title
#'
#' @param surv_fit_obj      survfit object 
#' @param time_unit         For the description in the table
#'
#' @returns                 A table with the time at 50% survival probability (+ confidnce intervalls) 
#'
#' @examples                Relevant for create_surv_plot()
get_surv_summary <- function(surv_fit_obj, time_unit = "Years") {
  
  # summary object
  tab <- summary(surv_fit_obj)$table
  
    # Get values
      records <- tab["records"]
      events  <- tab["events"]
      median  <- round(tab["median"] / 12, 2)
      lcl     <- round(tab["0.95LCL"] / 12, 2)
      ucl     <- round(tab["0.95UCL"] / 12, 2)
    
    # Output
    data.frame(
      #Records = records,
      #Events  = events,
      #Median  = median,
      #lcl     = lcl,
      #ucl     = ucl,
      Median_survival = paste0(median, "  [", lcl, "–", ucl, "]"),
      check.names = FALSE  # 👈 verhindert Umwandlung des Spaltennamens
    ) %>%
      setNames(paste0("Median\n Survival, ", time_unit,"\n[95%-CI]"))
}














# create table for survival values at specific times 
# Follow up Survival probability (in %)
#' Title
#'
#' @param surv_fit_obj     survfit object 
#' @param times            Default = 12, 24, 60 months
#' @param time_unit        For the description in the table
#'
#' @returns                A table with follow up survival probabilities (+ confidnce intervalls) 
#'                         at 1, 2, and 5 years
#'                         
#' @examples               relevant for create_surv_plot() 
get_surv_time <- function(surv_fit_obj, times = c(12, 24, 60), time_unit = "Years") {
  # summary object
  surv_summary <- summary(surv_fit_obj, times = times)
  
  time_months <- surv_summary$time
  surv_percent <- round(surv_summary$surv * 100, 2)
  ci_lower <- round(surv_summary$lower * 100, 2)
  ci_upper <- round(surv_summary$upper * 100, 2)
  
  
  # table
  tibble(
    FU_in_Years = time_months / 12,
    Survival =  paste0(surv_percent, "% [", ci_lower, "–", ci_upper, "]")) %>%
    
    rename_with(~ c(paste0("Follow up\nin ",time_unit), 
                    paste0("Survival probability (%),\n[95%-CI]")))
}














# create a Kaplan-Meier Curve
#' Title
#'
#' @param data                     data
#' @param surv_fit_obj             survfit object
#' @param group_factor_formula     Default = 1; Change to calculate for group variable 
#' @param xscale                   Default = month_year; from get_xscale(); adjustable
#' @param scale_end                To define x-achsis end point
#' @param scale_break              Default = 24 months; Break points 
#' @param title                    Plot title
#' @param title_size               Font size of the title
#' @param text_size                
#' @param risk_table_size          Font size of the risk table
#' @param tbl2_x                   x-coordinate in the plot for the get_surv_summary() table
#' @param time_unit                For the description in the table
#' @param show_tables              Decide whether the tables should be visible in the plot
#'
#' @returns
#' @export
#'
#' @examples
create_surv_plot <- function(data,
                             surv_fit_obj,
                             group_factor_formula = "1",
                             xscale = "m_y", 
                             scale_end, 
                             scale_break = 24,
                             title,
                             title_size = 11, 
                             text_size = 7,
                             risk_table_size = 2.7,
                             tbl2_x,
                             time_unit = "Years",
                             show_tables = TRUE) {
  
  # Get the calc for the x-scale
  xtrans <- get_xscale(xscale)
  
  breaksX <- seq(0, scale_end, by = scale_break)
  breaksX2 <- seq(0, scale_end, by = scale_break)
  breaksY <- c(0, 25, 50, 75, 100)
  

  
  plot_obj <- ggsurvplot(surv_fit_obj,
                         data = data,
                         risk.table = "nrisk_cumevents", 
                         surv.median.line = "hv",
                         ggtheme = theme_pubr(),
                         palette = "black", 
                         censor.shape = "ǀ",
                         conf.int = TRUE,
                         xlab = "Years",
                         ylab = "Survival probability (in %)",
                         title = title,
                         risk.table.y.text = TRUE,
                         break.time.by = scale_break,
                         risk.table.fontsize = risk_table_size
  )
  

  
  plot_obj$plot <- plot_obj$plot + 
    scale_x_continuous(breaks = breaksX,
                       labels = round(breaksX * xtrans, 2), expand = c(0.05, 0.1)) + 
    scale_y_continuous(labels = breaksY) +
    theme(axis.text.x = element_text(size = title_size),
          axis.title.x = element_text(size = title_size, face = "bold"),
          axis.text.y = element_text(size = title_size, face = "bold"),
          axis.title.y = element_text(size = title_size, face = "bold"),
          title = element_text(size = title_size, face = "bold")
    ) + theme(legend.position = "none")
  
  
  if (show_tables){
    tbl1 <- get_surv_time(surv_fit_obj, time_unit = time_unit)     # get the Followup Survival Probability in %
    tbl1 <- as.data.frame(tbl1)
    tbl2 <- get_surv_summary(surv_fit_obj, time_unit = time_unit)   # get the Median time
    tbl2 <- as.data.frame(tbl2)
    
    
    plot_obj$plot <- plot_obj$plot + 
      ggpp::annotate(geom = "table", x = 0,  y = 0, label = tbl1) + 
      ggpp::annotate(geom = "table", x = tbl2_x,  y = 0, label = tbl2)
  }
  
  
  plot_obj$table <- plot_obj$table +
    scale_x_continuous(breaks = breaksX2, labels = round(breaksX2 * xtrans, 2), expand = c(0.05, 0.1)) +
    xlab("") +
    theme(axis.text.x = element_text(size = title_size),
          axis.title.y = element_text(size = title_size),
          title = element_text(size = title_size, face = "bold"))
  

  return(plot_obj)
}





