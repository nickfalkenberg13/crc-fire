# crc-fire
Analysis of colorectal cancer (CRC) incidence and survival in relation to air pollution and wildfires.

## Overview of Scripts
**read_freq_tables.R**

Processes data files obtained from SEER queries for a given site in the large intestine (left colon, right colon, rectum, etc.) and creates separate tables for incidence and survival analysis.

Inputs:
- The site of the cancer case
  - From this, the script locates all data files corresponding to that site, reads them, and combines them.
  - These files are not currently on the Github repository.

Outputs:
- The script itself does not output anything. It just contains functions.
- The main output of the functions is a list containing two tables: one for incidence analysis and one for survival analysis

**add_zero_rows.R**

When downloading data from SEER, rows with zero case counts were excluded. These rows are needed for the incidence analysis since a group with no cases is still informative. This script adds in any rows containing combinations of covariates that were excluded in the data query process. All added rows are given a `count` value of 0.

**write_surv_inc_tables.R**

Uses the functions in `read_freq_tables.R` and `add_zero_rows.R` to create and save incidence and survival tables for a given site.

Inputs:
- The site of the cancer case, which is used as input for the functions in `read_freq_tables.R`.
- A vector of SEER registry regions. For this analysis, the registries were split into three regions:
  - Eastern and Southern (ES)
  - Pacific (Pac)
  - Southwest & Midwest (SWMW)
 
Outputs:
- Incidence table corresponding to the specified site.
- Survival table corresponding to the specified site.

**read_offset_tables.R**

Reads and processes data downloaded from CDC WONDER that contains offsets for each combination of covariates used in the incidence analysis.

Input: Offset files downloaded from CDC WONDER (not currently on this repository).

Output: A processed table containing offset values for combinations of covariates used in this analysis.

**read_exp_cov_files.R**

Reads and combines data files for exposures (wildfire and air pollutants) and covariates (various county-level metrics).

Inputs:
- Air pollution data sets
  - PM2.5
  - NO2
  - Ozone
- Wildfire data
- Time-varying county variables
- County-level moving metrics
- CDC Places county metrics
- SEER CRC county metrics
- All input files are under data/exp_cov_inputs

Output: A merged dataset containing all of the above variables matched to their respective counties.

**prepare_analytic_datasets.R**

Functions for combining exposure, covariate, and offset files with incidence tables and exposure and covariate files with survival tables.

**condense_offset_table.R**

Processes the offset table output from `read_offset_tables.R` to match the categorical groupings used in the SEER data. The output is a recoded offset table.
