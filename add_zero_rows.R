#source("read_freq_tables.R") # to read county codes

#data_inc <- fread("data/incidence_table-LEFT.csv")

make_grid <- function(rg){
  # keys for county fips
  county_codes_path <- paste("~/SEER Files/county_codes_", rg, ".txt", sep = "")
  county_codes <- read_county_codes(county_codes_path)
  
  # table with every combination of variable values
  grd <- expand.grid(code = county_codes$code,
                     year_rx = 0:21,
                     age_group = 0:3,
                     sex = 0:1,
                     marital_status = 0:2,
                     race = 0:5,
                     stage = c("distant", "insitu", "localized", "regional"))
  
  grd <- as.data.table(grd)
  
  # join county_fips
  setkey(grd, code)
  setkey(county_codes, code)
  grd <- grd[county_codes,]
  
  grd <- grd %>% 
    mutate(county_fips = as.numeric(county_fips))
  
  return(grd)
}

# grids <- lapply(c("SWMW", "Pac"), make_grid)
# grids <- rbindlist(grids)
# grids$count <- 0 # add count column
# grids$site <- "LEFT" # add site column
# 
# ## add the missing zero rows to incidence data
# key_cols <- setdiff(names(grids), c("count", "county_name"))
# 
# # perform anti-join
# zero_rows <- grids[!data_inc, on = key_cols]
# 
# # add zero rows in
# data_inc_full <- rbind(data_inc, zero_rows)


add_zero_rows <- function(data_inc, site, reg_list){
  grids <- lapply(reg_list, make_grid)
  grids <- rbindlist(grids)
  grids$count <- 0 # add count column
  grids$site <- site # add site column
  
  # making county_fips column in data_inc numeric
  # so it can be joined with grids
  data_inc <- data_inc %>% 
    mutate(county_fips = as.numeric(county_fips))
  
  ## add the missing zero rows to incidence data
  key_cols <- setdiff(names(grids), c("count", "county_name"))
  
  # perform anti-join
  zero_rows <- grids[!data_inc, on = key_cols]
  
  # add zero rows in
  data_inc_full <- rbind(data_inc, zero_rows)
  
  return(data_inc_full)
}