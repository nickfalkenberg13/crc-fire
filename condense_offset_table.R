#' Condense Offset Table
#' The purpose of this script is to process the offset
#' table that has been combined from CDC Wonder downloads.
#' 
#' The groupings/categories used in CDC Wonder may not match
#' those used in the outcome data downloaded from SEER.
#' 
#' To prepare the data for analysis, the offset table will
#' be manipulated so that covariate combinations will match those
#' in the SEER incidence tables. The result will be a population count
#' for every row in the incidence table.

library(dplyr)

offset <- read.csv("data/offset_table.csv")

### create new variables with codings that match incidence tables
# codes based on SEER query dictionaries
age_code <- stack(list("0" = c("< 1 year",
                               "1-4 years",
                               "5-9 years",
                               "10-14 years",
                               "15-19 years",
                               "20-24 years",
                               "25-29 years"),
                       "1" = c("30-34 years",
                               "35-39 years",
                               "40-44 years",
                               "45-49 years"),
                       "2" = c("50-54 years",
                               "55-59 years",
                               "60-64 years ",
                               "65-69 years",
                               "70-74 years"),
                       "3" = c("75-79 years",
                               "80-84 years",
                               "85+ years")))

race_code <- data.frame(values = c("White", "Black or African American",
                                   "American Indian or Alaska Native",
                                   "Asian or Pacific Islander"),
                        race_RECD = 0:3)

# adding newly coded columns either by joining or mutating
offset_new <- offset %>%
  left_join(age_code, by = c("age" = "values")) %>%
  rename(age_RECD = ind) %>% 
  left_join(race_code, by = c("race" = "values")) %>% 
  mutate(year_RECD = year - 2000,
         race_RECD = ifelse(ethnicity == "Hispanic or Latino", 4, race_RECD),
         sex_RECD = ifelse(sex == "Female", 1, 0))

# condensing by summing count for each identical set of age, race, year, and sex
offset_new <- offset_new %>% 
  aggregate(count ~ county_code + age_RECD + race_RECD + year_RECD + sex_RECD, FUN = "sum") %>% 
  rename(count_offset = count)

write.csv(offset_new, "data/offset_table_RECD.csv", row.names = F)
