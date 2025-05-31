# Load libraries
library(phytools)
library(fst)
library(data.table)
library(dplyr)
library(tidyverse)
library(magrittr)
library(MuMIn)
library(lubridate)
require(tidyr)
library(readr)
library(arrow)
library(sjstats)
library(corrplot)
library(readxl)
library(readr)
library(lme4)
library(tigris)
library(ggplot2)
library(sf)
setwd("/n/dominici_nsaph_l3/Lab/projects/adrd_hup_tspm/analysis")

load("../data/scratch/ADRD_hospitalization_2016-2018_clin.RData") #ADRD hospitalizations
load("../data/county_corr_to_national_clin_groups.RData") # county similarity to national HUP

setwd("/n/dominici_nsaph_l3/Lab/projects/adrd_hup_tspm/code")
load("../data/county_corr_to_national_clin_groups.RData")
shapefile_path <- "/n/dominici_nsaph_l3/Lab/projects/adrd_hup_tspm/data/input/cb_2018_us_county_500k/cb_2018_us_county_500k.shp"

county_data <- st_read(shapefile_path) %>% 
  st_transform(crs = 5070)%>%
  dplyr::select(c("STATEFP", "GEOID", "geometry")) %>%
  dplyr::rename(county = GEOID)

# Join shapefile
county_sk_sh = county_sk %>%
  left_join(county_data) %>% 
  filter(!STATEFP %in% c( "02","15", "66", "72", "60", "69", "78"))

# Add state
states_sf <- st_read("../data/input/tl_2018_us_state/tl_2018_us_state.shp") %>%
  dplyr::select(c("STUSPS", "STATEFP")) %>%
  st_drop_geometry()

county_sk_sh = county_sk_sh %>%
  left_join(states_sf)

# 1-MAD
county_sk_sh = county_sk_sh %>%
  mutate(mad_compl = 1 - mad)


# Prepare variables -------------------------------------------------------
# Load denom files
files_denom <- c("/n/dominici_nsaph_l3/Lab/projects/analytic/mbsf_medpar_denom/mbsf_medpar_denom_2016.parquet",
                 "/n/dominici_nsaph_l3/Lab/projects/analytic/mbsf_medpar_denom/mbsf_medpar_denom_2017.parquet",
                 "/n/dominici_nsaph_l3/Lab/projects/analytic/mbsf_medpar_denom/mbsf_medpar_denom_2018.parquet")
denom <- as.data.frame(rbindlist(lapply(files_denom,
                                        read_parquet)))%>%
  select(bene_id, zip, year, county)

fips_df = denom %>%
  select(zip, county, year)%>%
  group_by(zip, year)%>% 
  filter(row_number() == 1) %>%
  ungroup()%>%
  na.omit()

# Merge county to patients
unique_patients = ADRD_hospitalization %>%
  group_by(bene_id)%>%
  filter(row_number() == 1) %>%
  select(bene_id, race_rti, dual, orec, zip, state, year) %>%
  left_join(fips_df)

# Ratios of interest per county within the ADRD cohort
ratios_county <- unique_patients %>%
  group_by(county) %>%
  summarise(
    adrd_count = n(), # total patients with adrd for county
    # Ratio for each level of race_rti
    ratio_white = sum(race_rti == "1")/ adrd_count,
    ratio_black = sum(race_rti == "2")/ adrd_count,
    ratio_asian = sum(race_rti == "4")/ adrd_count,
    ratio_hispanic = sum(race_rti == "5")/ adrd_count,
    
    # Ratio dual == 1
    ratio_dual = sum(dual == 1)/ adrd_count,
    
    # Ratio disability
    ratio_disability = sum(orec %in% c(1, 3))/ adrd_count) %>%
  select(county, starts_with("ratio"), adrd_count)


# Ratios of interest at county level

# Residents in urban areas
nhgis_urban_spanish <- read_csv("/n/dominici_nsaph_l3/Lab/projects/adrd_hup_tspm/data/input/nhgis0003_ds258_2020_county.csv") %>%
  select(GEOCODE, U7I001, U7I002, U7I003, U7K001, U7K003) %>%
  mutate(urban = U7I002/ U7I001)%>% 
  mutate(rural = U7I003/ U7I001)%>%
  rename(county = GEOCODE)

# Residents with less than high school degree
education_data <- read_excel("../data/input/Education.xlsx", sheet = 1)[,c(1,52)]
colnames(education_data) = c("county", "perc_less_hs")
education_data$perc_less_hs = education_data$perc_less_hs/100

# Residents unemployed
unemployment_data <- read_excel("../data/input/Unemployment.xlsx", sheet = 1) %>%
  select(FIPS_Code, Unemployment_rate_2018)%>%
  rename(county = FIPS_Code) %>%
  mutate(Unemployment_rate_2018 = Unemployment_rate_2018/100)


# Ratio of adrd patients over total number of medicare patients
hosp_files <- c("/n/dominici_nsaph_l3/Lab/projects/analytic/mbsf_medpar_denom/medpar_hospitalizations_2016.parquet", 
                "/n/dominici_nsaph_l3/Lab/projects/analytic/mbsf_medpar_denom/medpar_hospitalizations_2017.parquet", 
                "/n/dominici_nsaph_l3/Lab/projects/analytic/mbsf_medpar_denom/medpar_hospitalizations_2018.parquet")


hosp = as.data.frame(rbindlist(lapply(hosp_files,
                                      read_parquet)))%>%
  filter(src_admsn_cd %in% c("7","1","2"))%>%
  select(bene_id, admission_date, discharge_date, diagnoses,adm_id)

colnames(hosp)

ratios_county_adrd = hosp %>%
  mutate(year = year(admission_date)) %>%
  left_join(denom)%>%
  group_by(bene_id)%>%
  filter(row_number() == 1) %>%
  ungroup() %>%
  group_by(county) %>%
  summarise(total_count = n())%>%
  left_join(ratios_county) %>%
  mutate(adrd_ratio = adrd_count/total_count)


# Create complete dataset 
reg_data = ratios_county_adrd %>%
  left_join(unemployment_data)%>%
  left_join(education_data) %>%
  left_join(nhgis_urban_spanish)

# save(reg_data, file= "../data/reg_data.RData")
# load("../data/reg_data.RData")

# subset on complete cases and merge with similarity
var = c("county", "total_count", "ratio_white", "ratio_black", "ratio_asian", "ratio_hispanic", "ratio_dual", "ratio_disability", "adrd_count", "adrd_ratio","Unemployment_rate_2018", "perc_less_hs", "urban", "rural")
sub = reg_data[complete.cases(reg_data[,var]),]
summary(sub)

sub = county_sk_sh %>%
  left_join(sub)

# Regression --------------------------------------------------------------

# Random skewers
lm_multivar_re =lmer(random_skewers ~ ratio_dual + Unemployment_rate_2018 + perc_less_hs +rural 
                     + ratio_black + ratio_hispanic + ratio_asian + adrd_ratio + ratio_disability +
                       (1|STUSPS), data = sub) 
summary(lm_multivar_re)
r.squaredGLMM(lm_multivar_re)
parameters::p_value(lm_multivar_re)

# 1 - MAD
lm_multivar_re =lmer(mad_compl ~ ratio_dual + Unemployment_rate_2018 + perc_less_hs +rural 
                     + ratio_black + ratio_hispanic + ratio_asian + adrd_ratio + ratio_disability +
                       (1|STUSPS), data = sub) 
summary(lm_multivar_re)
r.squaredGLMM(lm_multivar_re)



