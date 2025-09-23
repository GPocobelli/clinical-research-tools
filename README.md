
# clinical-research-tools

<br>

## Overview  

 * **Survival Analysis**:  
   Core functions to support time-to-event analyses (OS, PFS, etc.).  
   Includes helpers for time calculation, axis scaling, and tables.  

 * **Kaplan–Meier plots**:  
   High-level wrapper around `ggsurvplot()` with options for  
   axis transformation, risk tables, median survival and follow-up times.  

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
fit  <- survfit(surv ~ 1, data = lung)

# Plot with tables
create_surv_plot(lung, fit, xscale = "m_y", time_unit = "Years")
```

<br><br>



<img width="652" height="469" alt="Rplot" src="https://github.com/user-attachments/assets/166d0f89-9bdd-4be0-b898-726910f7e69e" />
