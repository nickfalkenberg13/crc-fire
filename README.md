# crc-fire
## Project Overview

Master's Thesis, Harvard TH Chan School of Public Health SM80 Epidemiology Program

**Title**: Ambient air pollution exposure and colorectal cancer incidence in wildfire-impacted areas in the United States

This project is a population-based cohort study of colorectal cancer (CRC) cases diagnosed between 2000 and 2020 throughout the United States. The aim was to examine associations between air pollution, wildfire exposure, and CRC incidence.

### Study Population

- Approx. 600,000 CRC cases diagnosed between 2000 and 2020.
- Cancer cases recorded in Surveillance, Epidemiology, and End Results (SEER) registries.
- Cases diagnosed in 609 counties and the state of Alaska.

### Exposures

Concentrations of the following air pollutants:

1. PM<sub>2.5</sub>
2. NO<sub>2</sub>
3. O<sub>3</sub>

Wildfire exposure, measured in three ways:

1. Number of wildfires intersecting a county.
2. Proportion of county area burned by wildfires.
3. Proportion of county population living in census tracts intersecting wildfires.

### Outcomes

Colorectal cancer, as defined using ICD codes C18 (excluding C18.1 for appendix), C19, and C20, diagnosed between 2000 and 2020.

### Data Sources

- SEER
  - CRC cases
  - Individual-level data of cases
  - County-level data of testing prevalence
  - Rural-Urban Continuum Codes
- Atmospheric Composition Analysis Group
  - PM<sub>2.5</sub>
- Monitoring Trends in Burn Severity
  - Wildfires
- PLACES
  - County-level data on health and social factors
- CDC Wonder
  - Population totals stratified by county, year, and other factors




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
  - PM<sub>2.5</sub>
  - NO<sub>2</sub>
  - O<sub>3</sub>
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

### Modeling scripts

**model_functions.R**

Contains helper functions to streamline model running and obtaining results.

**models_incidence_x.R**

Runs Poisson models for air pollutant *x*. May also contain code running exploratory data analysis.

**models_incidence_pm25_strat_x.R**

Runs models to assess three-way interaction between PM<sub>2.5</sub>, wildfires, and covariate *x*. First runs models stratified by *x* to get individual effect estimates for each level of *x*, then runs models to assess interaction.
