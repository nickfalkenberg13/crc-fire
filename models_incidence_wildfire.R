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

ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = numwf_3cat, y = count)) +
  geom_boxplot() + 
  labs(title = "Number of Wildfires and Case Counts")

ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = numwf_3cat, y = count)) +
  geom_violin() + 
  ylim(0, 20) +
  labs(title = "Number of Wildfires and Case Counts")

ggplot(data = d_full[sample(nrow(d_full), 5e4),] %>% filter(count > 0), aes(x = numwf_3cat, y = count)) +
  geom_boxplot() + 
  labs(title = "Number of Wildfires and Case Counts (Positive Case Counts Only")

ggplot(data = d_full[sample(nrow(d_full), 5e4),] %>% filter(count > 0), aes(x = numwf_4cat, y = count)) +
  geom_boxplot() + 
  labs(title = "Number of Wildfires and Case Counts (Positive Case Counts Only")

ggplot(data = d_full[sample(nrow(d_full), 5e4),] %>% filter(count > 0), aes(x = areawf_3cat, y = count)) +
  geom_boxplot() + 
  labs(title = "Area of Wildfires and Case Counts (Positive Case Counts Only")

ggplot(data = d_full[sample(nrow(d_full), 5e4),] %>% filter(count > 0), aes(x = areawf_4cat, y = count)) +
  geom_boxplot() + 
  labs(title = "Area of Wildfires and Case Counts (Positive Case Counts Only")

# same plot but as a violin
ggplot(data = d_full[sample(nrow(d_full), 5e4),] %>% filter(count > 0), aes(x = areawf_4cat, y = count)) +
  geom_violin() + 
  labs(title = "Area of Wildfires and Case Counts (Positive Case Counts Only")

ggplot(data = d_full[sample(nrow(d_full), 5e4),] %>% filter(count > 0), aes(x = popwf_3cat, y = count)) +
  geom_boxplot() + 
  labs(title = "Population Exposed to Wildfires and Case Counts (Positive Case Counts Only")

ggplot(data = d_full[sample(nrow(d_full), 5e4),] %>% filter(count > 0), aes(x = popwf_4cat, y = count)) +
  geom_boxplot() + 
  labs(title = "Population Exposed to Wildfires and Case Counts (Positive Case Counts Only")

## correlation between wildfire metrics

ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = numwf_3cat, y = areawf_3cat)) +
  geom_jitter() + 
  labs(title = "Number of Wildfires and Area of Wildfires")

ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = numwf_3cat, y = popwf_3cat)) +
  geom_jitter() + 
  labs(title = "Number of Wildfires and Population in Wildfires")

ggplot(data = d_full[sample(nrow(d_full), 5e4),], aes(x = areawf_3cat, y = popwf_3cat)) +
  geom_jitter() + 
  labs(title = "Area of Wildfires and Population in Wildfires")

# Modeling ----------------------------------------------------------------
exposures <- c("numwf_3cat", "numwf_4cat", "areawf_3cat", "areawf_4cat", "popwf_3cat", "popwf_4cat")

### basic models ###
covars_basic <- c("SEER.registry",
                  "age_RECD",
                  "race_RECD",
                  "sex_RECD",
                  "marital_status",
                  "year_RECD",
                  "offset(log(count_offset))")

results <- tibble()
for (exposure in exposures) {
  est <- get_effect_estimate(A = exposure, covars = covars_basic, d = d_full)
  results <- rbind(results, est)
}
write.csv(results,
          file = "model_results/wf_basic.csv",
          row.names = F)

### full models ###
covars_full <- c("SEER.registry",
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

results <- tibble()
for (exposure in exposures) {
  est <- get_effect_estimate(A = exposure, covars = covars_full, d = d_full)
  results <- rbind(results, est)
}
write.csv(results,
          file = "model_results/wf_full.csv",
          row.names = F)

#### Stratified by Race ####
# first run individual models for each race
race_list <- c("White",
               "Black",
               "American Indian/Alaska Native",
               "Asian/Pacific Islander",
               "Hispanic/Latino")

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

results <- tibble()
for (r in 0:4) {
  race <- race_list[r+1]
  print(paste("Race", race, sep = " = "))
  d_race <- d_full %>% filter(race_RECD == r)
  for (exposure in exposures) {
    est <- get_effect_estimate(A = exposure, covars = covars_no_race, d = d_race)
    est$race <- race
    results <- rbind(results, est)
  }
  print("----------------------------------------------------------")
}
write.csv(results,
          file = "model_results/wf_int_race.csv",
          row.names = F)

#### Assessing effect modification ####
for (exposure in exposures) {
  print(paste("Evaluating interaction for exposure", exposure))
  ## run model with interaction term
  int_term <- paste(exposure, "race_RECD", sep = "*") # set interaction term
  formula_int <- mod_formula(c(int_term, covars_no_race))
  mod_int <- glm(formula_int,
                 data = d_full,
                 family = poisson)
  
  ## run model without interaction
  formula_no_int <- mod_formula(c(exposure, covars_full))
  mod_no_int <- glm(formula_no_int,
                    data = d_full,
                    family = poisson)
  
  ## run likelihood ratio test
  print(lrtest(mod_no_int, mod_int))
  print("_____________________________________________________________")
}
