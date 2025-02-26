
# Libraries & Source Code ---------------------------------------------------------------

library(MASS) # for negative binomial modeling
library(pscl) # testing over-dispersion and running zero-inflated models
library(car) # for VIF
library(ggplot2)
library(lmtest) # for likelihood ratio testing
library(sandwich) # for robust standard errors
library(performance) # for checking zero inflation
library(vcdExtra) # for checking zero inflation

source("prepare_analytic_datasets.R")
source("model_functions.R")

# Data Prep ---------------------------------------------------------------

inc_LEFT <- fread("data/incidence_table-LEFT.csv")
inc_RIGHT <- fread("data/incidence_table-RIGHT.csv")
inc_RECTUM <- fread("data/incidence_table-RECTUM.csv")
inc_OTHER <- fread("data/incidence_table-OTHER.csv")

## create dataset
# combine files for each stage
d <- rbindlist(list(inc_LEFT, inc_RIGHT, inc_RECTUM, inc_OTHER))

# collapse by these columns; will not be used for now
sum_by_col <- setdiff(names(d), c("site", "count"))
d <- d[, .(count = sum(count)), by = sum_by_col]

# run data through process pipeline
d_full <- prepare_inc_data(d)
d_full %>% filter(count_offset == 0) %>% nrow() # no rows with offset 0

rm(inc_LEFT, inc_OTHER, inc_RECTUM, inc_RIGHT)

# PM2.5 IQR ---------------------------------------------------------------

exposure <- "pm25_iqr"
wildfire_vars <- c("popwf_3cat", "popwf_4cat", "numwf_3cat", "numwf_4cat", "areawf_3cat", "areawf_4cat")
strat_val_list <- c("insitu", "localized", "regional", "distant") # specifying order of stage values
# covariates, excluding stratifying variable
covars_no_strat <- c("SEER.registry",
                     "age_RECD",
                     "race_RECD",
                     "sex_RECD",
                     "marital_status",
                     "year_RECD",
                     "medhouseinc",
                     "DIABETES",
                     "obesity",
                     "LPA",
                     "CSMOKING",
                     "BINGE",
                     "urban",
                     "offset(log(count_offset))")

results <- data.frame()
for (r in strat_val_list) {
  print(paste("Stage", r, sep = " = "))
  d_strat <- d_full %>% filter(stage == r)
  for (wf in wildfire_vars) {
    d_wf_split <- split(d_strat, d_strat[[wf]])
    for (d_wf in d_wf_split){
      val <- sort(unique(d_wf[[wf]]))
      print(paste("Results for", wf, "=", val))
      est <- get_effect_estimate(A = exposure, covars = covars_no_strat, d = d_wf)
      est$stage <- r
      est$wf_subset_var <- wf
      est$wf_subset_val <- val
      results <- rbind(results, est)
    }
    print("-----------")
  }
  print("------------------------------------------------------------")
}
write.csv(results,
          file = paste("model_results/wf_int_", exposure, "_stage.csv", sep = ""),
          row.names = F)

#### Assessing effect modification ####
# (within a level of stage) is there EM between pm2.5 and wf?
for (r in strat_val_list) {
  print(paste("Stage", r, sep = " = "))
  for (wf in wildfire_vars) {
    print(paste("Evaluating interaction for", wf))
    ## run model with interaction term
    int_term <- paste(exposure, wf, sep = "*") # set interaction term
    formula_int <- mod_formula(c(int_term, covars_no_strat))
    mod_int <- glm(formula_int,
                   data = d_full,
                   family = poisson,
                   subset = (stage == r))
    
    ## run model without interaction
    formula_no_int <- mod_formula(c(exposure, wf, covars_no_strat))
    mod_no_int <- glm(formula_no_int,
                      data = d_full,
                      family = poisson,
                      subset = (stage == r))
    
    ## run likelihood ratio test
    print(lrtest(mod_no_int, mod_int))
    print("________________________________")
  }
  print("------------------------------------------------------------")
  print("------------------------------------------------------------")
}

strat_var <- "stage"

for (wf in wildfire_vars) {
  print(paste("Evaluating interaction for", wf))
  
  ## run model with three-way interaction term
  int_term <- paste(exposure, wf, strat_var, sep = "*") # set interaction term
  formula_int <- mod_formula(c(int_term, covars_no_strat))
  mod_int <- glm(formula_int,
                 data = d_full,
                 family = poisson)
  
  ## run model without three-way interaction
  int_two_ways <- c(paste(exposure, wf, sep = "*"),
                    paste(exposure, strat_var, sep = "*"),
                    paste(wf, strat_var, sep = "*"))
  formula_no_int <- mod_formula(c(int_two_ways, covars_no_strat))
  mod_no_int <- glm(formula_no_int,
                    data = d_full,
                    family = poisson)
  
  ## run likelihood ratio test
  print(lrtest(mod_no_int, mod_int))
  print("________________________________")
}
