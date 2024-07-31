source("read_freq_tables.R")
source("add_zero_rows.R")

site <- "RECTUM"
regs <- c("ES", "Pac", "SWMW")
# read and process files for given site
freq <- process_freq_data(site)

# create survival and incidence tables
freq_surv <- create_surv_table(freq_table = freq$survival)
freq_inc <- create_inc_table(freq_table = freq$incidence)

# add zero rows to incidence tables
freq_inc <- add_zero_rows(freq_inc, site, regs)

# set directory and file names for outputs
dir_surv <- paste("data/survival_table-", site, ".csv", sep = "")
dir_inc <- paste("data/incidence_table-", site, ".csv", sep = "")

# save outputs
fwrite(freq_surv, dir_surv)
fwrite(freq_inc, dir_inc)
