wonder_data_dir <- "~/Data/wonder_data/"
setwd(wonder_data_dir)

library(dplyr)

# read and process an individual file
process_data_file <- function(dat_file_name){
  offset_dat <- read.table(dat_file_name,
                           header = T)
  # remove excess columns
  drop_columns <- c("Yearly.July.1st.Estimates.Code", "Age.Group.Code", "Ethnicity.Code", "Race.Code")
  offset_dat <- offset_dat %>% select(!drop_columns)
  
  # handle Alaska special case, which does not have county column
  if(grepl("AK;", dat_file_name)){
    offset_dat <- offset_dat %>% 
      mutate(county = "Alaska",
             county_code = "02900",
             .before = "Yearly.July.1st.Estimates")
  }
  
  # rename remaining columns
  new_col_names <- c("county", "county_code", "year", "age", "ethnicity", "race", "count")
  names(offset_dat) <- new_col_names
  
  return(offset_dat)
}

# read every file and merge them into one
combine_offset_data <- function(sex = c("Male", "Female")){
  #sex <- ifelse(male, "Male", "Female")
  file_list <- list.files(pattern = paste("DATA ONLY.*", sex, sep = ""))
  
  # reading each file and accumulating the combined table
  offset_table <- data.frame()
  for (f in file_list) {
    d <- process_data_file(f)
    offset_table <- rbind(offset_table, d)
  }
  offset_table$sex <- sex # create the 'sex' variable
  
  return(offset_table)
}

# combine the male and female tables
offset_table_m <- combine_offset_data(sex = "Male")
offset_table_f <- combine_offset_data(sex = "Female")
offset_table <- rbind(offset_table_m, offset_table_f)

# save the results
setwd("~/crc-fire/data")
write.csv(offset_table, "offset_table.csv", row.names = F)
