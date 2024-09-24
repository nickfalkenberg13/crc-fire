
source("prepare_analytic_datasets.R")

library(MASS) # for negative binomial modeling
library(pscl) # testing over-dispersion and running zero-inflated models
library(car) # for VIF
library(ggplot2)
library(lmtest) # for likelihood ratio testing
library(sandwich) # for robust standard errors
library(performance) # for checking zero inflation
library(vcdExtra) # for checking zero inflation


# Data Prep ---------------------------------------------------------------

inc_LEFT <- fread("data/incidence_table-LEFT.csv")
inc_RIGHT <- fread("data/incidence_table-RIGHT.csv")
inc_RECTUM <- fread("data/incidence_table-RECTUM.csv")
inc_OTHER <- fread("data/incidence_table-RECTUM.csv")

# function for printing model results with robust SEs
model_sum_rb <- function(pm){
  # pm is a poisson model fit using glm
  model_tib <- broom::tidy(coeftest(pm, vcov. = vcovHC(pm, type = 'HC3')), 
                           conf.int = T)
  model_tib$estimate_exp <- exp(model_tib$estimate)
  model_tib$lower_exp <- exp(model_tib$conf.low)
  model_tib$upper_exp <- exp(model_tib$conf.high)
  print(model_tib[, c("term",
                      "estimate", 
                      "conf.low", 
                      "conf.high", 
                      "estimate_exp", 
                      "lower_exp", 
                      "upper_exp", 
                      "p.value")], n = nrow(model_tib), width = Inf)
}

# function for creating model formulas
mod_formula <- function(trms, resp_var = "count"){
  str_trms <- paste(trms, collapse = " + ")
  str_form <- paste(resp_var, str_trms, sep = " ~ ")
  
  return(formula(str_form))
}

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

ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = pm25_per10, y = count)) +
  geom_point() + 
  labs(title = "PM2.5 and Case Counts")

ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = age_RECD, y = count)) +
  geom_boxplot() + 
  labs(title = "Case Counts by Age Group")

ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = count)) +
  geom_histogram(binwidth = 1) + 
  labs(title = "Distribution of Case Counts")

ggplot(data = d_full, aes(x = year, y = pm25_per10)) +
  stat_summary(geom = "col", fun = median) + 
  labs(title = "Median PM2.5 by Year")

ggplot(data = d_full, aes(x = year, y = count)) +
  stat_summary(geom = "col", fun = sum) + 
  labs(title = "Total Cases by Year")

ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = pm25_per10, y = no2_per10)) +
  geom_point() +
  labs(title = "PM2.5 vs NO2")

ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = pm25_per10, y = ozone_per10)) +
  geom_point() +
  labs(title = "PM2.5 vs Ozone")

# PM2.5 and cases by wildfire category
ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = pm25_per10, y = count, colour = factor(numwf_3cat_noakhi))) +
  geom_point() +
  scale_color_manual(labels = c("No wildfires", "Below median", "Above median"), 
                     values = c("0" = "forestgreen", "1" = "goldenrod", "2" = "firebrick2")) +
  labs(color = "Wildfire Category", y = "Cases", x = "PM2.5 (per 10)")

ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = pm25_per10, y = count, colour = factor(numwf_4cat_noakhi))) +
  geom_point() +
  scale_color_manual(labels = c("No wildfires", "1st tertile", "2nd tertile", "3rd tertile"),
                     values = c("0" = "forestgreen", "1" = "goldenrod1", "2" = "peru", "3" = "firebrick2")) +
  labs(color = "Wildfire Category", y = "Cases", x = "PM2.5 (per 10)")

ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = pm25_per10, y = count/count_offset, colour = factor(numwf_4cat_noakhi))) +
  geom_point() +
  scale_color_manual(labels = c("No wildfires", "1st tertile", "2nd tertile", "3rd tertile"),
                     values = c("0" = "forestgreen", "1" = "goldenrod1", "2" = "peru", "3" = "firebrick2")) +
  labs(color = "Wildfire Category", y = "Rate", x = "PM2.5 (per 10)")

ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = pm25_per10, y = count/count_offset, colour = factor(numwf_5cat_noakhi))) +
  geom_point() +
  scale_color_manual(labels = c("No wildfires", "1st quartile", "2nd quartile", "3rd quartile", "4th quartile"),
                     values = c("0" = "forestgreen", 
                                "1" = "goldenrod1", 
                                "2" = "peru", 
                                "3" = "firebrick1", 
                                "4" = "firebrick4")) +
  labs(color = "Wildfire Category", y = "Rate", x = "PM2.5 (per 10)", title = "CRC Rate and PM2.5 Levels by Number of Wildfires")

# PM2.5 and cases by racial category
ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = race_RECD, y = count)) +
  geom_violin() +
  ylim(1, 20) +
  scale_x_discrete(labels = c("White", "Black", "American Indian/Alaska Native", "Asian", "Hispanic")) +
  labs(y = "Cases", x = "Race")

# PM2.5 Per 10 ------------------------------------------------------------

#### Basic model ####
# creating model formula
mod_terms_basic <- c("pm25_per10",
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
mod_terms_full <- c("pm25_per10",
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
##### population in wildfires #####
## 3-level variable
# 0 pop in wildfires
psn_per10_3cat_0 <- glm(formula_full, 
                        data = d_full, 
                        subset = which(d_full$popwf_3cat == 0),
                        family = poisson)
model_sum_rb(psn_per10_3cat_0)

# Below median
psn_per10_3cat_1 <- glm(formula_full, 
                        data = d_full, 
                        subset = which(d_full$popwf_3cat == 1),
                        family = poisson)
model_sum_rb(psn_per10_3cat_1)

# Above or equal to median
psn_per10_3cat_2 <- glm(formula_full, 
                        data = d_full, 
                        subset = which(d_full$popwf_3cat == 2),
                        family = poisson)
model_sum_rb(psn_per10_3cat_2)

# interaction model
mod_terms_int <- c("pm25_per10*popwf_3cat",
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
formula_int <- mod_formula(mod_terms_int)
psn_per10_3cat_int <- glm(formula_int,
                          data = d_full,
                          family = poisson)
model_sum_rb(psn_per10_3cat_int)

# test for interaction (using likelihood ratio test)
formula_no_int <- mod_formula(c("pm25_per10", 
                                "popwf_3cat",
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
                                "offset(log(count_offset))"))
mod_no_int <- glm(formula_no_int,
                  data = d_full,
                  family = poisson)
lrtest(mod_no_int, psn_per10_3cat_int) # fail to reject model without interaction
remove(mod_no_int, formula_no_int)

## 4-level variables
# 0 pop in wildfires
psn_per10_4cat_0 <- glm(formula_full, 
                        data = d_full, 
                        subset = (popwf_4cat == 0),
                        family = poisson)
model_sum_rb(psn_per10_4cat_0)

# below first tertile
psn_per10_4cat_1 <- glm(formula_full, 
                        data = d_full, 
                        subset = (popwf_4cat == 1),
                        family = poisson)
model_sum_rb(psn_per10_4cat_1)

# greater than or equal to first tertile and below 2nd tertile
psn_per10_4cat_2 <- glm(formula_full, 
                        data = d_full, 
                        subset = (popwf_4cat == 2),
                        family = poisson)
model_sum_rb(psn_per10_4cat_2)

# greater than or equal to second tertile
psn_per10_4cat_3 <- glm(formula_full, 
                        data = d_full, 
                        subset = (popwf_4cat == 3),
                        family = poisson)
model_sum_rb(psn_per10_4cat_3) %>% filter(term == "pm25_per10")

# interaction model
mod_terms_int <- c("pm25_per10*popwf_4cat",
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
formula_int <- mod_formula(mod_terms_int)
psn_per10_4cat_int <- glm(formula_int,
                          data = d_full,
                          family = poisson)

# test for interaction (using likelihood ratio test)
formula_no_int <- mod_formula(c("pm25_per10", 
                                "popwf_4cat",
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
                                "offset(log(count_offset))"))
mod_no_int <- glm(formula_no_int,
                  data = d_full,
                  family = poisson)
lrtest(mod_no_int, psn_per10_4cat_int) # fail to reject model without interaction
remove(mod_no_int, formula_no_int)

##### number of wildfires #####
## 3-level variable
# 0 wildfires
psn_per10_3cat_num_0 <- glm(formula_full,
                            data = d_full,
                            subset = (numwf_3cat == 0),
                            family = poisson)
model_sum_rb(psn_per10_3cat_num_0) %>% filter(term == "pm25_per10")

# Below median
psn_per10_3cat_num_1 <- glm(formula_full,
                            data = d_full,
                            subset = (numwf_3cat == 1),
                            family = poisson)
model_sum_rb(psn_per10_3cat_num_1) %>% filter(term == "pm25_per10")

# Above or equal to median
psn_per10_3cat_num_2 <- glm(formula_full,
                            data = d_full,
                            subset = (numwf_3cat == 2),
                            family = poisson)
model_sum_rb(psn_per10_3cat_num_2) %>% filter(term == "pm25_per10")

# interaction model
mod_terms_int <- c("pm25_per10*numwf_3cat",
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
formula_int <- mod_formula(mod_terms_int)
psn_per10_3cat_num_int <- glm(formula_int,
                              data = d_full,
                              family = poisson)

# test for interaction (using likelihood ratio test)
formula_no_int <- mod_formula(c("pm25_per10", 
                                "numwf_3cat",
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
                                "offset(log(count_offset))"))
mod_no_int <- glm(formula_no_int,
                  data = d_full,
                  family = poisson)
lrtest(mod_no_int, psn_per10_3cat_num_int) # fail to reject model without interaction
remove(mod_no_int, formula_no_int)

## 4-level variables
# 0 wildfires
psn_per10_4cat_num_0 <- glm(formula_full,
                            data = d_full,
                            subset = (numwf_4cat == 0),
                            family = poisson)
model_sum_rb(psn_per10_4cat_num_0) %>% filter(term == "pm25_per10")

# below first tertile
psn_per10_4cat_num_1 <- glm(formula_full,
                            data = d_full,
                            subset = (numwf_4cat == 1),
                            family = poisson)
model_sum_rb(psn_per10_4cat_num_1) %>% filter(term == "pm25_per10")

# greater than or equal to first tertile and below 2nd tertile
psn_per10_4cat_num_2 <- glm(formula_full,
                            data = d_full,
                            subset = (numwf_4cat == 2),
                            family = poisson)
model_sum_rb(psn_per10_4cat_num_2) %>% filter(term == "pm25_per10")

# greater than or equal to second tertile
psn_per10_4cat_num_3 <- glm(formula_full,
                            data = d_full,
                            subset = (numwf_4cat == 3),
                            family = poisson)
model_sum_rb(psn_per10_4cat_num_3) %>% filter(term == "pm25_per10")

# interaction model
mod_terms_int <- c("pm25_per10*numwf_4cat",
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
formula_int <- mod_formula(mod_terms_int)
psn_per10_4cat_num_int <- glm(formula_int,
                              data = d_full,
                              family = poisson)

# test for interaction (using likelihood ratio test)
formula_no_int <- mod_formula(c("pm25_per10", 
                                "numwf_4cat",
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
                                "offset(log(count_offset))"))
mod_no_int <- glm(formula_no_int,
                  data = d_full,
                  family = poisson)
lrtest(mod_no_int, psn_per10_4cat_num_int) # fail to reject model without interaction
remove(mod_no_int, formula_no_int)

##### area of wildfires #####
## 3-level variable
# 0 wildfires
psn_per10_3cat_area_0 <- glm(formula_full,
                             data = d_full,
                             subset = (areawf_3cat == 0),
                             family = poisson)
model_sum_rb(psn_per10_3cat_area_0) %>% filter(term == "pm25_per10")

# Below median
psn_per10_3cat_area_1 <- glm(formula_full,
                             data = d_full,
                             subset = (areawf_3cat == 1),
                             family = poisson)
model_sum_rb(psn_per10_3cat_area_1) %>% filter(term == "pm25_per10")

# Above or equal to median
psn_per10_3cat_area_2 <- glm(formula_full,
                             data = d_full,
                             subset = (areawf_3cat == 2),
                             family = poisson)
model_sum_rb(psn_per10_3cat_area_2) %>% filter(term == "pm25_per10")

# interaction model
mod_terms_int <- c("pm25_per10*areawf_3cat",
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
formula_int <- mod_formula(mod_terms_int)
psn_per10_3cat_area_int <- glm(formula_int,
                               data = d_full,
                               family = poisson)

# test for interaction (using likelihood ratio test)
formula_no_int <- mod_formula(c("pm25_per10", 
                                "areawf_3cat",
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
                                "offset(log(count_offset))"))
mod_no_int <- glm(formula_no_int,
                  data = d_full,
                  family = poisson)
lrtest(mod_no_int, psn_per10_3cat_area_int) # fail to reject model without interaction
remove(mod_no_int, formula_no_int)

## 4-level variables
# 0 wildfires
psn_per10_4cat_area_0 <- glm(formula_full,
                             data = d_full,
                             subset = (areawf_4cat == 0),
                             family = poisson)
model_sum_rb(psn_per10_4cat_area_0) %>% filter(term == "pm25_per10")

# below first tertile
psn_per10_4cat_area_1 <- glm(formula_full,
                             data = d_full,
                             subset = (areawf_4cat == 1),
                             family = poisson)
model_sum_rb(psn_per10_4cat_area_1) %>% filter(term == "pm25_per10")

# greater than or equal to first tertile and below 2nd tertile
psn_per10_4cat_area_2 <- glm(formula_full,
                             data = d_full,
                             subset = (areawf_4cat == 2),
                             family = poisson)
model_sum_rb(psn_per10_4cat_area_2) %>% filter(term == "pm25_per10")

# greater than or equal to second tertile
psn_per10_4cat_area_3 <- glm(formula_full,
                             data = d_full,
                             subset = (areawf_4cat == 3),
                             family = poisson)
model_sum_rb(psn_per10_4cat_area_3) %>% filter(term == "pm25_per10")

# interaction model
mod_terms_int <- c("pm25_per10*areawf_4cat",
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
formula_int <- mod_formula(mod_terms_int)
psn_per10_4cat_area_int <- glm(formula_int,
                               data = d_full,
                               family = poisson)

# test for interaction (using likelihood ratio test)
formula_no_int <- mod_formula(c("pm25_per10", 
                                "areawf_4cat",
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
                                "offset(log(count_offset))"))
mod_no_int <- glm(formula_no_int,
                  data = d_full,
                  family = poisson)
lrtest(mod_no_int, psn_per10_4cat_area_int)
remove(mod_no_int, formula_no_int)


# PM2.5 IQR ---------------------------------------------------------------

#### Basic model ####
# creating model formula
exposure <- "pm25_iqr"
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
##### population in wildfires #####
## 3-level variable
# 0 pop in wildfires
psn_iqr_3cat_0 <- glm(formula_full,
                      data = d_full,
                      subset = (d_full$popwf_3cat == 0),
                      family = poisson)
model_sum_rb(psn_iqr_3cat_0) %>% filter(term == exposure)

# Below median
psn_iqr_3cat_1 <- glm(formula_full,
                      data = d_full,
                      subset = (d_full$popwf_3cat == 1),
                      family = poisson)
model_sum_rb(psn_iqr_3cat_1) %>% filter(term == exposure)

# Above or equal to median
psn_iqr_3cat_2 <- glm(formula_full,
                      data = d_full,
                      subset = (d_full$popwf_3cat == 2),
                      family = poisson)
model_sum_rb(psn_iqr_3cat_2) %>% filter(term == exposure)

# interaction model
wildfire_var <- "popwf_3cat"
mod_terms_int <- c(paste(exposure, wildfire_var, sep = "*"),
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
formula_int <- mod_formula(mod_terms_int)
psn_iqr_3cat_int <- glm(formula_int,
                        data = d_full,
                        family = poisson)

# test for interaction (using likelihood ratio test)
formula_no_int <- mod_formula(c(exposure, 
                                wildfire_var,
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
                                "offset(log(count_offset))"))
mod_no_int <- glm(formula_no_int,
                  data = d_full,
                  family = poisson)
lrtest(mod_no_int, psn_iqr_3cat_int) # fail to reject model without interaction
remove(mod_no_int, formula_no_int)

## 4-level variables
# 0 pop in wildfires
psn_iqr_4cat_0 <- glm(formula_full,
                      data = d_full,
                      subset = (popwf_4cat == 0),
                      family = poisson)
model_sum_rb(psn_iqr_4cat_0) %>% filter(term == exposure)

# below first tertile
psn_iqr_4cat_1 <- glm(formula_full,
                      data = d_full,
                      subset = (popwf_4cat == 1),
                      family = poisson)
model_sum_rb(psn_iqr_4cat_1) %>% filter(term == exposure)

# greater than or equal to first tertile and below 2nd tertile
psn_iqr_4cat_2 <- glm(formula_full,
                      data = d_full,
                      subset = (popwf_4cat == 2),
                      family = poisson)
model_sum_rb(psn_iqr_4cat_2) %>% filter(term == exposure)

# greater than or equal to second tertile
psn_iqr_4cat_3 <- glm(formula_full,
                      data = d_full,
                      subset = (popwf_4cat == 3),
                      family = poisson)
model_sum_rb(psn_iqr_4cat_3) %>% filter(term == exposure)

# interaction model
wildfire_var <- "popwf_4cat"
mod_terms_int <- c(paste(exposure, wildfire_var, sep = "*"),
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
formula_int <- mod_formula(mod_terms_int)
psn_iqr_4cat_int <- glm(formula_int,
                        data = d_full,
                        family = poisson)

# test for interaction (using likelihood ratio test)
formula_no_int <- mod_formula(c(exposure, 
                                wildfire_var,
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
                                "offset(log(count_offset))"))
mod_no_int <- glm(formula_no_int,
                  data = d_full,
                  family = poisson)
lrtest(mod_no_int, psn_iqr_4cat_int)
remove(mod_no_int, formula_no_int)

##### number of wildfires #####
## 3-level variable
# 0 wildfires
psn_iqr_3cat_num_0 <- glm(formula_full,
                      data = d_full,
                      subset = (d_full$numwf_3cat == 0),
                      family = poisson)
model_sum_rb(psn_iqr_3cat_num_0) %>% filter(term == exposure)

# Below median
psn_iqr_3cat_num_1 <- glm(formula_full,
                      data = d_full,
                      subset = (d_full$numwf_3cat == 1),
                      family = poisson)
model_sum_rb(psn_iqr_3cat_num_1) %>% filter(term == exposure)

# Above or equal to median
psn_iqr_3cat_num_2 <- glm(formula_full,
                      data = d_full,
                      subset = (d_full$numwf_3cat == 2),
                      family = poisson)
model_sum_rb(psn_iqr_3cat_num_2) %>% filter(term == exposure)

# interaction model
wildfire_var <- "numwf_3cat"
mod_terms_int <- c(paste(exposure, wildfire_var, sep = "*"),
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
formula_int <- mod_formula(mod_terms_int)
psn_iqr_3cat_num_int <- glm(formula_int,
                        data = d_full,
                        family = poisson)

# test for interaction (using likelihood ratio test)
formula_no_int <- mod_formula(c(exposure, 
                                wildfire_var,
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
                                "offset(log(count_offset))"))
mod_no_int <- glm(formula_no_int,
                  data = d_full,
                  family = poisson)
lrtest(mod_no_int, psn_iqr_3cat_num_int) # fail to reject model without interaction
remove(mod_no_int, formula_no_int)

## 4-level variables
# 0 wildfires
psn_iqr_4cat_num_0 <- glm(formula_full,
                          data = d_full,
                          subset = (numwf_4cat == 0),
                          family = poisson)
model_sum_rb(psn_iqr_4cat_num_0) %>% filter(term == exposure)

# below first tertile
psn_iqr_4cat_num_1 <- glm(formula_full,
                          data = d_full,
                          subset = (numwf_4cat == 1),
                          family = poisson)
model_sum_rb(psn_iqr_4cat_num_1) %>% filter(term == exposure)

# greater than or equal to first tertile and below 2nd tertile
psn_iqr_4cat_num_2 <- glm(formula_full,
                          data = d_full,
                          subset = (numwf_4cat == 2),
                          family = poisson)
model_sum_rb(psn_iqr_4cat_num_2) %>% filter(term == exposure)

# greater than or equal to second tertile
psn_iqr_4cat_num_3 <- glm(formula_full,
                          data = d_full,
                          subset = (numwf_4cat == 3),
                          family = poisson)
model_sum_rb(psn_iqr_4cat_num_3) %>% filter(term == exposure)

# interaction model
wildfire_var <- "numwf_4cat"
mod_terms_int <- c(paste(exposure, wildfire_var, sep = "*"),
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
formula_int <- mod_formula(mod_terms_int)
psn_iqr_4cat_num_int <- glm(formula_int,
                            data = d_full,
                            family = poisson)

# test for interaction (using likelihood ratio test)
formula_no_int <- mod_formula(c(exposure, 
                                wildfire_var,
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
                                "offset(log(count_offset))"))
mod_no_int <- glm(formula_no_int,
                  data = d_full,
                  family = poisson)
lrtest(mod_no_int, psn_iqr_4cat_num_int) # fail to reject model without interaction
remove(mod_no_int, formula_no_int)

##### area of wildfires #####
## 3-level variable
# 0 wildfires
psn_iqr_3cat_area_0 <- glm(formula_full,
                          data = d_full,
                          subset = (d_full$areawf_3cat == 0),
                          family = poisson)
model_sum_rb(psn_iqr_3cat_area_0) %>% filter(term == exposure)

# Below median
psn_iqr_3cat_area_1 <- glm(formula_full,
                          data = d_full,
                          subset = (d_full$areawf_3cat == 1),
                          family = poisson)
model_sum_rb(psn_iqr_3cat_area_1) %>% filter(term == exposure)

# Above or equal to median
psn_iqr_3cat_area_2 <- glm(formula_full,
                          data = d_full,
                          subset = (d_full$areawf_3cat == 2),
                          family = poisson)
model_sum_rb(psn_iqr_3cat_area_2) %>% filter(term == exposure)

# interaction model
wildfire_var <- "areawf_3cat"
mod_terms_int <- c(paste(exposure, wildfire_var, sep = "*"),
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
formula_int <- mod_formula(mod_terms_int)
psn_iqr_3cat_area_int <- glm(formula_int,
                            data = d_full,
                            family = poisson)

# test for interaction (using likelihood ratio test)
formula_no_int <- mod_formula(c(exposure, 
                                wildfire_var,
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
                                "offset(log(count_offset))"))
mod_no_int <- glm(formula_no_int,
                  data = d_full,
                  family = poisson)
lrtest(mod_no_int, psn_iqr_3cat_area_int) # fail to reject model without interaction
remove(mod_no_int, formula_no_int)

## 4-level variables
# 0 wildfires
psn_iqr_4cat_area_0 <- glm(formula_full,
                          data = d_full,
                          subset = (areawf_4cat == 0),
                          family = poisson)
model_sum_rb(psn_iqr_4cat_area_0) %>% filter(term == exposure)

# below first tertile
psn_iqr_4cat_area_1 <- glm(formula_full,
                          data = d_full,
                          subset = (areawf_4cat == 1),
                          family = poisson)
model_sum_rb(psn_iqr_4cat_area_1) %>% filter(term == exposure)

# greater than or equal to first tertile and below 2nd tertile
psn_iqr_4cat_area_2 <- glm(formula_full,
                          data = d_full,
                          subset = (areawf_4cat == 2),
                          family = poisson)
model_sum_rb(psn_iqr_4cat_area_2) %>% filter(term == exposure)

# greater than or equal to second tertile
psn_iqr_4cat_area_3 <- glm(formula_full,
                          data = d_full,
                          subset = (areawf_4cat == 3),
                          family = poisson)
model_sum_rb(psn_iqr_4cat_area_3) %>% filter(term == exposure)

# interaction model
wildfire_var <- "areawf_4cat"
mod_terms_int <- c(paste(exposure, wildfire_var, sep = "*"),
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
formula_int <- mod_formula(mod_terms_int)
psn_iqr_4cat_area_int <- glm(formula_int,
                            data = d_full,
                            family = poisson)

# test for interaction (using likelihood ratio test)
formula_no_int <- mod_formula(c(exposure, 
                                wildfire_var,
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
                                "offset(log(count_offset))"))
mod_no_int <- glm(formula_no_int,
                  data = d_full,
                  family = poisson)
lrtest(mod_no_int, psn_iqr_4cat_area_int) # fail to reject model without interaction
remove(mod_no_int, formula_no_int)
