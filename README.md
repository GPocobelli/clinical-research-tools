
# clinical-research-tools

<br>

## Overview  

 * **Survival Analysis**:  
   Core functions to support time-to-event analyses (OS, PFS, etc.).  
   Includes helpers for time calculation, axis scaling, and tables.  

 * **Kaplan–Meier plots**:  
   High-level wrapper around `ggsurvplot()` with options for  
   axis transformation, risk tables, median survival and follow-up times.

   > :warning: Just for one group for now

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
    x_end = 1000,                # optional: end of the x-axis in days
    title = "Kaplan–Meier Curve",
    tbl2_pos = c(0.8, 0.99),     # position of the follow up table
    risk_table_size = 3.5        # fontsize of the values
)
```

<br><br>

![Rplot](https://github.com/user-attachments/assets/eab49722-653b-4ec2-b758-baa821b97d41)
