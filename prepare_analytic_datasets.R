#' Prepare Analytic Datasets
#' 
#' This script prepares incidence and survival data for modeling by
#' combining it with covariate/exposure data and offset
#' data from CDC Wonder. It then applies other necessary
#' data manipulations and returns the output.
#' 
#' The functions in this script are meant to be used after
#' outcome data is subsetted as necessary (e.g. looking at
#' specific sites or stages in the stratified analysis).
#' 

library(dplyr)
library(data.table)

#### Incidence Data ####
# # incidence data
# data_inc <- fread("data/incidence_table-LEFT.csv")
# # exposure and covariate data
# data_exp_co <- fread("data/merged_exposure_covariate_data.csv")
# # offset data
# data_offset <- fread("data/offset_table_RECD.csv")
# 
# setkey(data_inc, county_fips, year_rx, age_group, sex, race)
# setkey(data_offset, county_code, year_RECD, age_RECD, sex_RECD, race_RECD)
# data_exp_co$year_RECD <- data_exp_co$year - 2000
# setkey(data_exp_co, countyfips, year_RECD)
# # join offset and incidence
# data_cmbn <- data_offset[data_inc, on = .(county_code = county_fips,
#                                           year_RECD = year_rx,
#                                           age_RECD = age_group,
#                                           sex_RECD = sex,
#                                           race_RECD = race)]
# 
# # add covariates
# data_cmbn <- data_exp_co[data_cmbn, on =.(countyfips = county_code,
#                                           year_RECD)]
# 
# # assigning NA values to 0 for offset
# # these are rows that did not show up in the offset table
# # since zero values were filtered in the CDC Wonder queries
# data_cmbn[is.na(data_cmbn$count_offset), "count_offset"] <- 0
# 
# fwrite(data_cmbn, "data/analytic_dataset_incidence.csv")

prepare_inc_data <- function(data_inc){
  # exposure and covariate data
  data_exp_co <- fread("data/merged_exposure_covariate_data.csv")
  # offset data
  data_offset <- fread("data/offset_table_RECD.csv")
  
  setkey(data_inc, county_fips, year_rx, age_group, sex, race)
  setkey(data_offset, county_code, year_RECD, age_RECD, sex_RECD, race_RECD)
  data_exp_co$year_RECD <- data_exp_co$year - 2000
  setkey(data_exp_co, countyfips, year_RECD)
  # join offset and incidence
  data_cmbn <- data_offset[data_inc, on = .(county_code = county_fips,
                                            year_RECD = year_rx,
                                            age_RECD = age_group,
                                            sex_RECD = sex,
                                            race_RECD = race)]
  
  # add covariates
  data_cmbn <- data_exp_co[data_cmbn, on =.(countyfips = county_code,
                                            year_RECD)]
  
  # assigning NA values to 0 for offset
  # these are rows that did not show up in the offset table
  # since zero values were filtered in the CDC Wonder queries
  data_cmbn[is.na(data_cmbn$count_offset), "count_offset"] <- 0
  
  ## removing rows with missing data and special cases
  # cases in county 35999 are duplicates of Cibola and Valencia Counties
  # unknown counties cannot have exposure and covariates assigned to them
  # 2021 is excluded in Wonder, so no offset would be available
  # race category 5 is unknown race
  data_cmbn <- data_cmbn %>%
    filter(countyfips != 35999 &
             !grepl("Unknown", county_name) &
             year_RECD < 21 &
             race_RECD != 5)
  
  return(data_cmbn)
}

#### Survival Data ####
prepare_surv_data <- function(data_surv){
  # exposure and covariate data
  data_exp_co <- fread("data/merged_exposure_covariate_data.csv")
  
  data_exp_co$year_RECD <- data_exp_co$year - 2000
  setkey(data_exp_co, countyfips, year_RECD)
  
  # add covariates
  data_cmbn <- data_exp_co[data_surv, on =.(countyfips = county_fips,
                                            year_RECD = year_rx)]
  
  return(data_cmbn)
}