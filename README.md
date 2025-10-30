
# clinical-research-tools

<br>

## Overview  

 * **Survival Analysis**:  
   Core functions to support time-to-event analyses (OS, PFS, etc.).  
   Includes helpers for time calculation, axis scaling, and tables.  

 * **Kaplan–Meier plots**:  
   High-level wrapper around `ggsurvplot()` with options for  
   axis transformation, risk tables, median survival and follow-up times.


<br>

 * **Tables**:  
   Ready-to-annotate tables for median survival and follow-up probabilities  
   (e.g. at 1, 2, 5 years).  

<br>
<br>  

 

## Procedure  

### 1. Ensure the following structure exists  
 - Package dependencies: **survival**, **survminer**, **tidyverse**, **ggpubr**, **ggpp**.  
 - Functions are available via the script **`survival_utils.R`**:  

```r
source("path/to/survival_util.R")
```

<br> <br>

## Available functions  

- calc_time() – difference between two dates in days, months, or years  
- get_xscale() – factor for transforming the time axis
- get_median_table() – one-row table of median survival with 95% CI
- get_surv_times() – survival probabilities at predefined follow-up times
- create_surv_plot() – wrapper to produce Kaplan–Meier plots with tables


<br><br>

## Typical usage
### One group

```r
# Kaplan–Meier fit
library(survival) 
surv <- Surv(time = lung$time, event = lung$status)
fit <- survfit(surv ~ 1, data = lung)

create_surv_plot(
    data = lung,
    fit = fit,
    xscale = "d_m",              # days to months
    time_unit = "Months",        
    scale_break = 182.625,       # 6 months (182.625 days)
    title = "Kaplan–Meier Curve",
    tbl1_pos = c(0, 1.2),
    tbl2_pos = c(1, 1.2),     # position of the follow up table
    risk_table_size = 3.5        # fontsize of the values
)
```

<br><br>

<img width="976" height="631" alt="image" src="https://github.com/user-attachments/assets/d5b0017f-640f-4b96-995c-8cf010f4d2c8" />



<br><br><br><br>



### Two groups

```r
# Kaplan–Meier fit
lung <- survival::lung %>%
    mutate(sex = case_when(sex == 1 ~ "Male", sex == 2 ~ "Female"))
fit <- survfit(Surv(time, status) ~ sex, data = lung)

create_surv_plot(lung,
                 fit,
                 xscale = "d_m",
                 scale_break = 182.625,
                 title = "Kaplan–Meier Curve",
                 tbl1_pos = c(0.3, 1.2),
                 tbl2_pos = c(1, 1.2),
                 risk.table.y.text = FALSE,
                 risk.table.y.text.col = TRUE,
                 risk.table.col = "black")
```

<br><br>

<img width="976" height="604" alt="image" src="https://github.com/user-attachments/assets/e2eb911d-07c3-4c05-843b-5e9f09a9be86" />

