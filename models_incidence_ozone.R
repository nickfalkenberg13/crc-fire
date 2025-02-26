
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

## Ozone and cases by wildfire category ##
# num_wf
ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = ozone_per10, y = count, colour = factor(numwf_3cat))) +
  geom_point() +
  scale_color_manual(labels = c("No wildfires", "Below median", "Above median"), 
                     values = c("0" = "forestgreen", "1" = "goldenrod", "2" = "firebrick2")) +
  labs(color = "Wildfire Category", y = "Cases", x = "Ozone (per 10)")

ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = ozone_iqr, y = count, colour = factor(numwf_3cat))) +
  geom_point() +
  scale_color_manual(labels = c("No wildfires", "Below median", "Above median"), 
                     values = c("0" = "forestgreen", "1" = "goldenrod", "2" = "firebrick2")) +
  labs(color = "Wildfire Category", y = "Cases", x = "Ozone (per IQR)")

# area_wf
ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = ozone_per10, y = count, colour = factor(areawf_3cat))) +
  geom_point() +
  scale_color_manual(labels = c("No wildfires", "Below median", "Above median"), 
                     values = c("0" = "forestgreen", "1" = "goldenrod", "2" = "firebrick2")) +
  labs(color = "Wildfire Category", y = "Cases", x = "Ozone (per 10)")

ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = ozone_iqr, y = count, colour = factor(areawf_3cat))) +
  geom_point() +
  scale_color_manual(labels = c("No wildfires", "Below median", "Above median"), 
                     values = c("0" = "forestgreen", "1" = "goldenrod", "2" = "firebrick2")) +
  labs(color = "Wildfire Category", y = "Cases", x = "Ozone (per IQR)")

## Ozone and case rate by wildfire category ##
# num_wf
ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = ozone_per10, y = count/count_offset, colour = factor(numwf_3cat))) +
  geom_point() +
  scale_color_manual(labels = c("No wildfires", "Below median", "Above median"), 
                     values = c("0" = "forestgreen", "1" = "goldenrod", "2" = "firebrick2")) +
  labs(color = "Number of Wildfires", y = "Rate", x = "Ozone (per 10)")

# area_wf
ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = ozone_per10, y = count/count_offset, colour = factor(areawf_3cat))) +
  geom_point() +
  scale_color_manual(labels = c("No wildfires", "Below median", "Above median"), 
                     values = c("0" = "forestgreen", "1" = "goldenrod", "2" = "firebrick2")) +
  labs(color = "Area of Wildfires", y = "Rate", x = "Ozone (per 10)")

# pop_wf
ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = ozone_per10, y = count/count_offset, colour = factor(popwf_3cat))) +
  geom_point() +
  scale_color_manual(labels = c("No wildfires", "Below median", "Above median"), 
                     values = c("0" = "forestgreen", "1" = "goldenrod", "2" = "firebrick2")) +
  labs(color = "Population in Wildfires", y = "Rate", x = "Ozone (per 10)")

## Ozone levels by wildfire category ##
# num_wf
ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = numwf_3cat, y = ozone_per10)) +
  geom_boxplot() + 
  labs(title = "Ozone levels by Wildfire Category", x = "Wildfire Category (# of Wildfires)", y = "Ozone (per 10)")

ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = numwf_3cat, y = ozone_iqr)) +
  geom_boxplot() + 
  labs(title = "Ozone levels by Wildfire Category", x = "Wildfire Category (# of Wildfires)", y = "Ozone (per IQR)")

# area_wf
ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = areawf_3cat, y = ozone_per10)) +
  geom_boxplot() + 
  labs(title = "Ozone levels by Wildfire Category", x = "Wildfire Category (Area of Wildfires)", y = "Ozone (per 10)")

ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = areawf_4cat, y = ozone_per10)) +
  geom_boxplot() + 
  labs(title = "Ozone levels by Wildfire Category", x = "Wildfire Category (Area of Wildfires)", y = "Ozone (per 10)")

ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = areawf_4cat, y = ozone_iqr)) +
  geom_boxplot() + 
  labs(title = "Ozone levels by Wildfire Category", x = "Wildfire Category (Area of Wildfires)", y = "Ozone (per IQR)")

## Cases by racial category ##
ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = race_RECD, y = count)) +
  geom_violin() +
  ylim(1, 20) +
  scale_x_discrete(labels = c("White", "Black", "American Indian/Alaska Native", "Asian", "Hispanic")) +
  labs(y = "Cases", x = "Race")

ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = race_RECD, y = count)) +
  geom_boxplot() +
  ylim(1, 20) +
  scale_x_discrete(labels = c("White", "Black", "American Indian/Alaska Native", "Asian", "Hispanic")) +
  labs(y = "Cases", x = "Race")


# Ozone per 10 ------------------------------------------------------------

exposure <- "ozone_per10"

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
exposure <- "ozone_per10"
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
  formula_int <- mod_formula(c(int_term, mod_terms_full[-1]))
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


# Ozone IQR -----------------------------------------------------------

exposure <- "ozone_iqr"

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
model_sum_rb(psn_iqr_basic)

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
model_sum_rb(psn_iqr_full)

#### Stratified by wildfire category ####

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
exposure <- "ozone_iqr"
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
  formula_int <- mod_formula(c(int_term, mod_terms_full[-1]))
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
