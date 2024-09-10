library(dplyr)

setwd("data/exp_cov_inputs")

#### Read files ####
### Exposures
# wildfire data file
data_exp_wf <- read.csv("wfvars_final.csv")

## Pollutants
# PM 2.5 data file 
data_exp_pm25 <- read.csv("pmvars_final.csv")

# NO2 data file 
data_exp_no2 <- read.csv("no2vars_final.csv")

# Ozone data file 
data_exp_ozn <- read.csv("ozonevars_final.csv")

### Covariates
# Time-varying county variables data file
data_cnty_tv <- read.csv("countyvars_tv_final.csv")

# Moving-related county variables data file
data_cnty_move <- read.csv("countyvars_moving_final.csv")

# CDC Places county data file
data_cnty_places <- read.csv("PLACES_FINAL.csv")

# SEER CRC County Metrics
data_cnty_crc <- read.csv("crc_county_metrics.csv")

#### Join data ####
### Time varying data (per county per year)
## NO2 and Ozone
data_cnty_yr <- full_join(data_exp_no2, data_exp_ozn, 
                          by = c("countyfips", "SEER.registry", "GEOID10", "year"))
## PM 2.5
# introduces rows for AK and HI, which will be NA for NO2 and Ozone columns
data_cnty_yr <- full_join(data_cnty_yr, data_exp_pm25,
                          by = c("countyfips", "SEER.registry", "GEOID10", "year"))

## Wildfires
data_cnty_yr <- full_join(data_cnty_yr, data_exp_wf,
                          by = c("countyfips", "SEER.registry", "GEOID10", "year" = "year_wf"))

## County variables
# doing a left join because the county variables have extra years and counties that are not needed
data_cnty_yr <- left_join(data_cnty_yr, data_cnty_tv,
                          by = c("countyfips", "year"))

### Time constant data (per county only)
data_cnty_only <- full_join(data_cnty_move, data_cnty_places, by = "countyfips")

## Preparing CRC metrics
names(data_cnty_crc) <- c("CountyName", 
                          "pct_CRC_2004-2007",
                          "pct_Endo_2004-2007",
                          "pct_FOBT_2004-2007",
                          "pct_CRC_2008-2010",
                          "pct_Endo_2008-2010",
                          "pct_FOBT_2008-2010")

# Extracting FIPS and converting metrics columns to numbers
data_cnty_crc <- data_cnty_crc %>% 
  filter(grepl("[0-9]+", CountyName)) %>% 
  mutate(countyfips = regmatches(CountyName,
                                  regexpr("[0-9]+", CountyName)),
         countyfips = as.numeric(countyfips),
         `pct_CRC_2004-2007` = as.numeric(`pct_CRC_2004-2007`)/100,
         `pct_Endo_2004-2007` = as.numeric(`pct_Endo_2004-2007`)/100,
         `pct_FOBT_2004-2007` = as.numeric(`pct_FOBT_2004-2007`)/100,
         `pct_CRC_2008-2010` = as.numeric(`pct_CRC_2008-2010`)/100,
         `pct_Endo_2008-2010` = as.numeric(`pct_Endo_2008-2010`)/100,
         `pct_FOBT_2008-2010` = as.numeric(`pct_FOBT_2008-2010`)/100)

data_cnty_only <- left_join(data_cnty_only, data_cnty_crc, by = "countyfips")

### Combine everything
data_merged <- left_join(data_cnty_yr, data_cnty_only, by = "countyfips", relationship = "many-to-one")

# write results
write.csv(data_merged, "../merged_exposure_covariate_data.csv", row.names = F)
