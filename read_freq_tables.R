library(dplyr)
library(data.table)

#site <- "LEFT"

initial_surv_mo <- function(ft){
  # returns the first number in the range of survival months covered in  frequency table ft
  
  surv_col <- names(ft)[grepl("Survival months", names(ft))] # column name containing survival months
  month_int <- strsplit(surv_col, split = " ")[[1]][3]
  first_month <- strsplit(month_int, split = "-")[[1]][1]
  
  return(as.numeric(first_month))
}

read_freq_file <- function(fpath, stg, ste){
  #cases <- read.csv(fpath)
  cases <- fread(fpath)
  #cases <- freq_file %>% filter(Count != 0)
  initial_m <- initial_surv_mo(cases)
  names(cases) <- c("county", "year_rx", "age_group", "sex", "marital_status",
                    "race", "cod", "survival_months", "count")
  
  # add initial survival month to survival month column to get accurate survival numbers
  # add stage and site as columns
  cases <- cases %>% mutate(survival_months = survival_months + initial_m,
                            stage = stg,
                            site = ste)
  
  return(cases)
}

read_county_codes <- function(fpath){
  cnty_codes <- read.table(fpath,
                           sep = "=",
                           col.names = c("code", "county_name"))
  
  # creating column of FIPS codes, extracted using REGEX
  cnty_codes <- cnty_codes %>% mutate(county_fips = regmatches(county_name,
                                                               regexpr("[0-9]+", county_name)))
  return(cnty_codes)
}

process_freq_data <- function(ste){
  table_path <- paste("~/SEER Files/frequency_site", ste, sep = "-")
  #stages <- c("distant", "insitu", "localized", "regional")
  #registry_groups <- c("ES", "Pac", "SWMW")
  stages <- c("distant", "insitu")
  registry_groups <- c("Pac", "SWMW")
  
  data_output <- data.frame()
  
  for (rg in registry_groups) {
    # read and combine files for all registry groups for the site ste
    data_registry <- data.frame()
    
    for (st in stages) {
      print(paste("Processing files of registry", rg, "and stage", st))
      
      # read and combine files of all stages for a registry group
      file_list <- list.files(path = table_path, 
                              pattern = paste("freq", rg, st, "csv", sep = ".*"))
      data_stage <- data.frame()
      
      for (f in file_list) {
        # read and combine all files of a given stage and registry group
        f_path <- paste(table_path, f, sep = "/")
        cases <- read_freq_file(f_path, st, ste)
        data_stage <- rbind(data_stage, cases)
        
        print(paste("Finished processing file", f))
      }
      
      data_registry <- rbind(data_registry, data_stage)
    }
    
    # get county codes for given registries
    county_codes_path <- paste("~/SEER Files/county_codes_", rg, ".txt", sep = "")
    county_codes <- read_county_codes(county_codes_path)
    
    # joining the FIPS to the frequency table
    data_registry <- data_registry %>% 
      left_join(county_codes, by = join_by(county == code))
    
    data_output <- rbind(data_output, data_registry)
    
    print(paste("Finished processing", rg, "files."))
    print("-----------------------")
    print("-----------------------")
  }
  
  # moving the count variable to the end
  data_output <- data_output %>% relocate(count, .after = last_col())
  
  return(data_output)
}

#freq_left <- process_freq_data(site)
#write.csv(freq_left, "data/frequency_table_left.csv", row.names = F)

create_surv_table <- function(freq_table){
  # removes rows with 0 counts
  return(freq_table %>% filter(count > 0))
}

create_inc_table <- function(freq_table){
  freq_dt <- as.data.table(freq_table)
  
  # removes survival_months and cause of death columns
  # and sums counts of otherwise identical rows
  sum_by_col <- setdiff(names(freq_dt), c("survival_months", "cod", "count"))
  result <- freq_dt[, .(count = sum(count)), by = sum_by_col]
  
  return(result)
}
