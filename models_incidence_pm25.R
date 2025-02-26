
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
lrtest(mod_no_int, psn_iqr_3cat_int)
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

#### doing all of the above in a loop ####
wildfire_vars <- c("popwf_3cat", "popwf_4cat", "numwf_3cat", "numwf_4cat", "areawf_3cat", "areawf_4cat")
exposure <- "pm25_iqr"
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

results <- data.frame()
for (wf in wildfire_vars) {
  d_wf_split <- split(d_full, d_full[[wf]])
  for (d_wf in d_wf_split){
    val <- sort(unique(d_wf[[wf]]))
    print(paste("Results for", wf, "=", val))
    est <- get_effect_estimate(A = exposure, covars = mod_terms_full[-1], d = d_wf)
    est$wf_subset_var <- wf
    est$wf_subset_val <- val
    results <- rbind(results, est)
  }
  print("-----------")
}
write.csv(results,
          file = paste("model_results/wf_int_", exposure, ".csv", sep = ""),
          row.names = F)

# assess interaction
for (wf in wildfire_vars) {
  print(paste("Evaluating interaction for", wf))
  ## run model with interaction term
  int_term <- paste(exposure, wf, sep = "*") # set interaction term
  formula_int <- mod_formula(c(int_term, mod_terms_full[-1]))
  mod_int <- glm(formula_int,
                 data = d_full,
                 family = poisson)
  
  ## run model without interaction
  formula_no_int <- mod_formula(c(wf, mod_terms_full))
  mod_no_int <- glm(formula_no_int,
                    data = d_full,
                    family = poisson)
  
  ## run likelihood ratio test
  print(lrtest(mod_no_int, mod_int))
  print("________________________________")
}

# Stratified by Race ------------------------------------------------------

# Repeat above wildfire stratified analysis for each race
wildfire_vars <- c("popwf_3cat", "popwf_4cat", "numwf_3cat", "numwf_4cat", "areawf_3cat", "areawf_4cat")
race_list <- c("White",
               "Black",
               "American Indian/Alaska Native",
               "Asian/Pacific Islander",
               "Hispanic/Latino")

## PM2.5 IQR ---------------------------------------------------------------
exposure <- "pm25_iqr"
covars_no_race <- c("SEER.registry",
                    "age_RECD",
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
for (r in 0:4) {
  race <- race_list[r+1]
  print(paste("Race", race, sep = " = "))
  d_race <- d_full %>% filter(race_RECD == r)
  for (wf in wildfire_vars) {
    d_wf_split <- split(d_race, d_race[[wf]])
    for (d_wf in d_wf_split){
      val <- sort(unique(d_wf[[wf]]))
      print(paste("Results for", wf, "=", val))
      est <- get_effect_estimate(A = exposure, covars = covars_no_race, d = d_wf)
      est$race <- race
      est$wf_subset_var <- wf
      est$wf_subset_val <- val
      results <- rbind(results, est)
    }
    print("-----------")
  }
  print("------------------------------------------------------------")
}
write.csv(results,
          file = paste("model_results/wf_int_", exposure, "_race.csv", sep = ""),
          row.names = F)

#### Assessing effect modification ####
# (within a level of race) is there EM between pm2.5 and wf?
for (r in 0:4) {
  print(paste("Race", race_list[r+1], sep = " = "))
  for (wf in wildfire_vars) {
    print(paste("Evaluating interaction for", wf))
    ## run model with interaction term
    int_term <- paste(exposure, wf, sep = "*") # set interaction term
    formula_int <- mod_formula(c(int_term, covars_no_race))
    mod_int <- glm(formula_int,
                   data = d_full,
                   family = poisson,
                   subset = (race_RECD == r))
    
    ## run model without interaction
    formula_no_int <- mod_formula(c(exposure, wf, covars_no_race))
    mod_no_int <- glm(formula_no_int,
                      data = d_full,
                      family = poisson,
                      subset = (race_RECD == r))
    
    ## run likelihood ratio test
    print(lrtest(mod_no_int, mod_int))
    print("________________________________")
  }
  print("------------------------------------------------------------")
  print("------------------------------------------------------------")
}

#### Three-way interaction ####
# comparing model with three way interaction (and sub-interactions)
# to model with no interactions
for (wf in wildfire_vars) {
  print(paste("Evaluating interaction for", wf))
  
  ## run model with three-way interaction term
  int_term <- paste(exposure, wf, "race_RECD", sep = "*") # set interaction term
  formula_int <- mod_formula(c(int_term, covars_no_race))
  mod_int <- glm(formula_int,
                 data = d_full,
                 family = poisson)
  
  ## run model without three-way interaction
  int_two_ways <- c(paste(exposure, wf, sep = "*"),
                    paste(exposure, "race_RECD", sep = "*"),
                    paste(wf, "race_RECD", sep = "*"))
  formula_no_int <- mod_formula(c(int_two_ways, covars_no_race))
  mod_no_int <- glm(formula_no_int,
                    data = d_full,
                    family = poisson)
  
  ## run likelihood ratio test
  print(lrtest(mod_no_int, mod_int))
  print("________________________________")
}
