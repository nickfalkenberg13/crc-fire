#' Prepare Analytic Datasets
#' 
#' This script prepares the final datasets that will
#' be used for modeling. It combines the merged covariate dataset
#' with the outcome data from SEER and the offset data from CDC Wonder.
#' 
#' It then applies any other remaining data manipulations required for
#' the analysis. Finally, it saves the output.
#' 

library(dplyr)
library(data.table)

#### Incidence Data ####
# incidence data
data_inc <- fread("data/incidence_table-LEFT.csv")
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

fwrite(data_cmbn, "data/analytic_dataset_incidence.csv")
