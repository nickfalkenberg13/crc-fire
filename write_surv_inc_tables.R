source("read_freq_tables.R")

site <- "LEFT"
# read and process files for given site
freq <- process_freq_data(site)

# create survival and incidence tables
freq_surv <- create_surv_table(freq_table = freq)
freq_inc <- create_inc_table(freq_table = freq)

# set directory and file names for outputs
dir_surv <- paste("data/survival_table-", site, ".csv", sep = "")
dir_inc <- paste("data/incidence_table-", site, ".csv", sep = "")

# save outputs
fwrite(freq_surv, dir_surv)
fwrite(freq_inc, dir_inc)
