
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

## PM2.5 and cases by wildfire category and sex ##
# num_wf
ggplot(data = d_full[sample(nrow(d_full), 5e4),], 
       aes(x = pm25_per10, y = count, colour = factor(numwf_3cat), shape = sex_RECD)) +
  geom_point() +
  scale_color_manual(labels = c("No wildfires", "Below median", "Above median"), 
                     values = c("0" = "forestgreen", "1" = "goldenrod", "2" = "firebrick2")) +
  scale_shape_manual(labels = c("Male", "Female"),
                     values = 0:1) +
  labs(color = "Number of Wildfires", y = "Cases", x = "PM2.5 (per 10)", shape = "Sex")

## PM2.5 and cases by sex ##
ggplot(data = d_full[sample(nrow(d_full), 5e4),], 
       aes(x = pm25_per10, y = count, color = sex_RECD)) +
  geom_point() +
  scale_color_manual(labels = c("Male", "Female"),
                     values = c("0" = "turquoise", "1" = "plum")) +
  labs(y = "Cases", x = "PM2.5 (per 10)", color = "Sex")

## Cases by sex ##
ggplot(data = d_full[sample(nrow(d_full), 5e4),], 
       aes(x = sex_RECD, y = count)) +
  geom_boxplot() +
  labs(y = "Cases", x = "Sex")


# PM2.5 per 10 ------------------------------------------------------------

exposure <- "pm25_per10"
wildfire_vars <- c("popwf_3cat", "popwf_4cat", "numwf_3cat", "numwf_4cat", "areawf_3cat", "areawf_4cat")
strat_val_list <- 0:1 # sex values
# covariates, excluding stratifying variable
covars_no_strat <- c("SEER.registry",
                     "age_RECD",
                     "race_RECD",
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
  sex <- c("Male", "Female")[r+1]
  print(paste("Sex", sex, sep = " = "))
  d_strat <- d_full %>% filter(sex_RECD == r)
  for (wf in wildfire_vars) {
    d_wf_split <- split(d_strat, d_strat[[wf]])
    for (d_wf in d_wf_split){
      val <- sort(unique(d_wf[[wf]]))
      print(paste("Results for", wf, "=", val))
      est <- get_effect_estimate(A = exposure, covars = covars_no_strat, d = d_wf)
      est$sex <- sex
      est$wf_subset_var <- wf
      est$wf_subset_val <- val
      results <- rbind(results, est)
    }
    print("-----------")
  }
  print("------------------------------------------------------------")
}
write.csv(results,
          file = paste("model_results/wf_int_", exposure, "_sex.csv", sep = ""),
          row.names = F)

#### Assessing effect modification ####
# (within a level of sex) is there EM between pm2.5 and wf?
for (r in strat_val_list) {
  print(paste("Sex", c("Male", "Female")[r+1], sep = " = "))
  for (wf in wildfire_vars) {
    print(paste("Evaluating interaction for", wf))
    ## run model with interaction term
    int_term <- paste(exposure, wf, sep = "*") # set interaction term
    formula_int <- mod_formula(c(int_term, covars_no_strat))
    mod_int <- glm(formula_int,
                   data = d_full,
                   family = poisson,
                   subset = (sex_RECD == r))
    
    ## run model without interaction
    formula_no_int <- mod_formula(c(exposure, wf, covars_no_strat))
    mod_no_int <- glm(formula_no_int,
                      data = d_full,
                      family = poisson,
                      subset = (sex_RECD == r))
    
    ## run likelihood ratio test
    print(lrtest(mod_no_int, mod_int))
    print("________________________________")
  }
  print("------------------------------------------------------------")
  print("------------------------------------------------------------")
}

#### Three-way interaction ####
# comparing model with three way interaction (and sub-interactions)
# to model with no three way interaction
strat_var <- "sex_RECD"

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


# PM2.5 IQR ---------------------------------------------------------------

exposure <- "pm25_iqr"
wildfire_vars <- c("popwf_3cat", "popwf_4cat", "numwf_3cat", "numwf_4cat", "areawf_3cat", "areawf_4cat")
strat_val_list <- 0:1 # sex values
# covariates, excluding stratifying variable
covars_no_strat <- c("SEER.registry",
                     "age_RECD",
                     "race_RECD",
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
  sex <- c("Male", "Female")[r+1]
  print(paste("Sex", sex, sep = " = "))
  d_strat <- d_full %>% filter(sex_RECD == r)
  for (wf in wildfire_vars) {
    d_wf_split <- split(d_strat, d_strat[[wf]])
    for (d_wf in d_wf_split){
      val <- sort(unique(d_wf[[wf]]))
      print(paste("Results for", wf, "=", val))
      est <- get_effect_estimate(A = exposure, covars = covars_no_strat, d = d_wf)
      est$sex <- sex
      est$wf_subset_var <- wf
      est$wf_subset_val <- val
      results <- rbind(results, est)
    }
    print("-----------")
  }
  print("------------------------------------------------------------")
}
write.csv(results,
          file = paste("model_results/wf_int_", exposure, "_sex.csv", sep = ""),
          row.names = F)

#### Assessing effect modification ####
# (within a level of sex) is there EM between pm2.5 and wf?
for (r in strat_val_list) {
  print(paste("Sex", c("Male", "Female")[r+1], sep = " = "))
  for (wf in wildfire_vars) {
    print(paste("Evaluating interaction for", wf))
    ## run model with interaction term
    int_term <- paste(exposure, wf, sep = "*") # set interaction term
    formula_int <- mod_formula(c(int_term, covars_no_strat))
    mod_int <- glm(formula_int,
                   data = d_full,
                   family = poisson,
                   subset = (sex_RECD == r))
    
    ## run model without interaction
    formula_no_int <- mod_formula(c(exposure, wf, covars_no_strat))
    mod_no_int <- glm(formula_no_int,
                      data = d_full,
                      family = poisson,
                      subset = (sex_RECD == r))
    
    ## run likelihood ratio test
    print(lrtest(mod_no_int, mod_int))
    print("________________________________")
  }
  print("------------------------------------------------------------")
  print("------------------------------------------------------------")
}

#### Three-way interaction ####
# comparing model with three way interaction (and sub-interactions)
# to model with no three way interaction
strat_var <- "sex_RECD"

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
