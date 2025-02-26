
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

# function for printing model results with robust SEs
# model_sum_rb <- function(pm){
#   # pm is a poisson model fit using glm
#   model_tib <- broom::tidy(coeftest(pm, vcov. = vcovHC(pm, type = 'HC3')), 
#                            conf.int = T)
#   model_tib$estimate_exp <- exp(model_tib$estimate)
#   model_tib$lower_exp <- exp(model_tib$conf.low)
#   model_tib$upper_exp <- exp(model_tib$conf.high)
#   # print(model_tib[, c("term",
#   #                     "estimate", 
#   #                     "conf.low", 
#   #                     "conf.high", 
#   #                     "estimate_exp", 
#   #                     "lower_exp", 
#   #                     "upper_exp", 
#   #                     "p.value")], n = nrow(model_tib), width = Inf)
#   return(model_tib[, c("term",
#                        "estimate",
#                        "conf.low",
#                        "conf.high",
#                        "estimate_exp",
#                        "lower_exp",
#                        "upper_exp",
#                        "p.value")])
# }
# 
# # function for creating model formulas
# mod_formula <- function(trms, resp_var = "count"){
#   str_trms <- paste(trms, collapse = " + ")
#   str_form <- paste(resp_var, str_trms, sep = " ~ ")
#   
#   return(formula(str_form))
# }

## create dataset
# combine files for each site
d <- rbindlist(list(inc_LEFT, inc_RIGHT, inc_RECTUM, inc_OTHER))
# collapse by these columns; will not be used for now
sum_by_col <- setdiff(names(d), c("stage", "site", "count"))
d <- d[, .(count = sum(count)), by = sum_by_col]

# run data through process pipeline
d_full <- prepare_inc_data(d)
d_full %>% filter(count_offset == 0) %>% nrow() # no rows with offset 0

rm(inc_LEFT, inc_OTHER, inc_RECTUM, inc_RIGHT)

# Exploratory Analysis ----------------------------------------------------


# NO2 Per 10 ------------------------------------------------------------
exposure <- "no2_per10"

#### Basic model ####
# creating model formula
mod_terms_basic <- c(exposure,
                     "SEER.registry",
                     "age_RECD",
                     "race_RECD",
                     "sex_RECD",
                     "marital_status",
                     "year_RECD",
                     "offset(log(count_offset))")
formula__basic <- mod_formula(mod_terms_basic)

# training model
psn_per10_basic <- glm(formula__basic, data = d_full, family = poisson)
model_sum_rb(psn_per10_basic)

#### Fully-adjusted model #####
# creating model formula
mod_terms_full <- c(exposure,
                    "SEER.registry",
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
formula_full <- mod_formula(mod_terms_full)

# training model
psn_per10_full <- glm(formula_full, data = d_full, family = poisson)
model_sum_rb(psn_per10_full)

#### Stratified by wildfire category ####
## will try to streamline this as done with stratified race analysis
## in models_incidence_pm25.R

# function to run models
# run_poisson <- function(A, covars, d, return_model = F){
#   mod_terms <- c(A, covars)
#   frmla <- mod_formula(mod_terms)
#   
#   psn_mod <- glm(frmla, data = d, family = poisson)
#   print(paste("Model results for exposure", A))
#   print(model_sum_rb(psn_mod) %>% filter(grepl(A, term)))
#   print("--------------------------")
#   
#   if(return_model){
#     return(psn_mod)
#   }
# }

wildfire_vars <- c("numwf_3cat", "numwf_4cat",
                   "areawf_3cat", "areawf_4cat",
                   "popwf_3cat", "popwf_4cat")

# run individual models to get effect estimates
results <- data.frame()
for (wf in wildfire_vars) {
  d_wf_split <- split(d_full, d_full[[wf]])
  for (d_wf in d_wf_split){
    val <- sort(unique(d_wf[[wf]]))
    print(paste("Results for", wf, "=", val))
    est <- get_effect_estimate(A = exposure, covars = mod_terms_full[-1], d = d_wf)
    params <- data.frame("Wf_Subset_Var" = wf,
                         "Wf_Subset_Val" = val)
    est <- cbind(est, params)
    results <- rbind(results, est)
  }
  print("-----------")
}
write.csv(results,
          file = paste("model_results/wf_int_", exposure, ".csv", sep = ""),
          row.names = F)

## test for interaction
# make sure we have the right exposure
exposure <- "no2_per10"
mod_terms_full <- c(exposure,
                    "SEER.registry",
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

# conduct likelihood ratio test to assess interaction
for (wf in wildfire_vars) {
  print(paste("Evaluating interaction for", wf))
  ## run model with interaction term
  int_term <- paste(exposure, wf, sep = "*") # set interaction term
  formula_int <- mod_formula(c(int_term, mod_terms_full))
  mod_int <- glm(formula_int,
                 data = d_full,
                 family = poisson)
  
  ## run model without interaction
  # exposure + wildfire as individual model terms
  formula_no_int <- mod_formula(c(wf, mod_terms_full))
  mod_no_int <- glm(formula_no_int,
                    data = d_full,
                    family = poisson)
  
  ## run likelihood ratio test
  print(lrtest(mod_no_int, mod_int))
  print("________________________________")
}

# NO2 IQR ---------------------------------------------------------------
exposure <- "no2_iqr"

#### Basic model ####
# creating model formula
mod_terms_basic <- c(exposure,
                     "SEER.registry",
                     "age_RECD",
                     "race_RECD",
                     "sex_RECD",
                     "marital_status",
                     "year_RECD",
                     "offset(log(count_offset))")
formula__basic <- mod_formula(mod_terms_basic)

# training model
psn_iqr_basic <- glm(formula__basic, data = d_full, family = poisson)
model_sum_rb(psn_iqr_basic) %>% filter(term == exposure)

#### Fully-adjusted model #####
# creating model formula
mod_terms_full <- c(exposure,
                    "SEER.registry",
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
formula_full <- mod_formula(mod_terms_full)

# training model
psn_iqr_full <- glm(formula_full, data = d_full, family = poisson)
model_sum_rb(psn_iqr_full) %>% filter(term == exposure)

#### Stratified by wildfire category ####
results <- data.frame()
for (wf in wildfire_vars) {
  d_wf_split <- split(d_full, d_full[[wf]])
  for (d_wf in d_wf_split){
    val <- sort(unique(d_wf[[wf]]))
    print(paste("Results for", wf, "=", val))
    est <- get_effect_estimate(A = exposure, covars = mod_terms_full[-1], d = d_wf)
    params <- data.frame("Wf_Subset_Var" = wf,
                         "Wf_Subset_Val" = val)
    est <- cbind(est, params)
    results <- rbind(results, est)
  }
  print("-----------")
}
write.csv(results,
          file = paste("model_results/wf_int_", exposure, ".csv", sep = ""),
          row.names = F)

## test for interaction
# make sure we have the right exposure
exposure <- "no2_iqr"
mod_terms_full <- c(exposure,
                    "SEER.registry",
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

# conduct likelihood ratio test to assess interaction
for (wf in wildfire_vars) {
  print(paste("Evaluating interaction for", wf))
  ## run model with interaction term
  int_term <- paste(exposure, wf, sep = "*") # set interaction term
  formula_int <- mod_formula(c(int_term, mod_terms_full))
  mod_int <- glm(formula_int,
                 data = d_full,
                 family = poisson)
  
  ## run model without interaction
  # exposure + wildfire as individual model terms
  formula_no_int <- mod_formula(c(wf, mod_terms_full))
  mod_no_int <- glm(formula_no_int,
                    data = d_full,
                    family = poisson)
  
  ## run likelihood ratio test
  print(lrtest(mod_no_int, mod_int))
  print("________________________________")
}
