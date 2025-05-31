# Required libraries
library(table1)
library(raster)
library(foreign)
library(fst)
library(data.table)
library(dplyr)
library(tidyverse)
library(magrittr)
library(lubridate)
library(readr)
library(arrow)
library(readxl)
setwd("/n/dominici_nsaph_l3/Lab/projects/adrd_hup_tspm/analysis")

# Load denom and hospitalizations files
files_denom <- c("/n/dominici_nsaph_l3/Lab/projects/analytic/mbsf_medpar_denom/mbsf_medpar_denom_2016.parquet",
                 "/n/dominici_nsaph_l3/Lab/projects/analytic/mbsf_medpar_denom/mbsf_medpar_denom_2017.parquet",
                 "/n/dominici_nsaph_l3/Lab/projects/analytic/mbsf_medpar_denom/mbsf_medpar_denom_2018.parquet")
denom <- as.data.frame(rbindlist(lapply(files_denom,
                                        read_parquet)))%>%
  select(bene_id, zip, year, race, race_rti, age_dob, state, dual, sex, orec)

hosp_files <- c("/n/dominici_nsaph_l3/Lab/projects/analytic/mbsf_medpar_denom/medpar_hospitalizations_2016.parquet", 
                "/n/dominici_nsaph_l3/Lab/projects/analytic/mbsf_medpar_denom/medpar_hospitalizations_2017.parquet", 
                "/n/dominici_nsaph_l3/Lab/projects/analytic/mbsf_medpar_denom/medpar_hospitalizations_2018.parquet")

# Filter on hospitalizations not coming from nursing facilities
hosp = as.data.frame(rbindlist(lapply(hosp_files,
                                      read_parquet)))%>%
  filter(src_admsn_cd %in% c("7","1","2"))%>%
  select(bene_id, admission_date, discharge_date, diagnoses,adm_id)

# Adjust data
hosp = hosp %>% unnest(diagnoses) %>% group_by(adm_id) %>% 
  mutate(col=seq_along(adm_id)) %>%
  spread(key=col, value=diagnoses)
  
colnames(hosp)[5:29] = paste0("DIAG", 1:25)

## Identify ADRD codes ----
# Change variables in character
codes_var = paste0("DIAG", 1:25)
hosp[codes_var] <- lapply(hosp[codes_var], as.character) 

# ADRD codes 
dementia_codes <- read.csv("../data/ADRD_ICD10.csv", stringsAsFactors = FALSE)
names(dementia_codes) <- c("ICD10", "category", "Description")
ADRD_codes = dementia_codes$ICD10

# Identify and select ADRD hospitalization: 
# select patients with at least one qualifying claim among the first 10 billing codes
ADRD_hospitalization = hosp %>% 
  filter(any(across(starts_with("DIAG"), ~`%in%`(.,ADRD_codes))))

# Join individual data
ADRD_hospitalization = ADRD_hospitalization %>%
  mutate(year = year(admission_date)) %>%
  left_join(denom, by =c ("bene_id", "year"))

# Filter on age < 100
ADRD_hospitalization = ADRD_hospitalization %>% 
  filter(age_dob < 100)

# Filter on patients that stayed in the same state during study 
# Same state 
n_states = ADRD_hospitalization %>%
  group_by(bene_id) %>%
  summarise(n_states = n_distinct(state))

bene_states = n_states[n_states$n_states==1, "bene_id"]
ADRD_hospitalization = ADRD_hospitalization[ADRD_hospitalization$bene_id %in% bene_states$bene_id, ]

# Same zipcode 
n_zipcodes = ADRD_hospitalization %>%
  group_by(bene_id) %>%
  summarise(n_zipcodes = n_distinct(zip))

bene_zip = n_zipcodes[n_zipcodes$n_zipcodes==1, "bene_id"]
ADRD_hospitalization = ADRD_hospitalization[ADRD_hospitalization$bene_id %in% bene_zip$bene_id, ]

save(ADRD_hospitalization, file = "../data/scratch/ADRD_hospitalization_2016-2018_clin.RData")


## Table 1 ---------------------------------------------------
load("../data/scratch/ADRD_hospitalization_2016-2018_clin.RData")

# Number of hospitalizations and patients
length(unique(ADRD_hospitalization$bene_id))
nrow(ADRD_hospitalization)

# Avg number of hospitalization per patient
mean(table(ADRD_hospitalization$bene_id))
sd(table(ADRD_hospitalization$bene_id))

# Table 1 
colnames(denom)
colnames(ADRD_hospitalization)

denom <- denom %>%
  select(bene_id, zip, year, race, race_rti, age_dob, state, dual, sex)

data_tab1 <- ADRD_hospitalization %>%
  ungroup() %>%
  select(bene_id, year) %>%
  group_by(bene_id, year)%>%
  filter(row_number() == 1) %>%
  left_join(denom)%>%
  group_by(bene_id)%>%
  filter(row_number() == 1)


summary(data_tab1)
table(data_tab1$sex)

mean(data_tab1$age_dob)
sd(data_tab1$age_dob)

#### CREATING TABLE 1######
data_tab1$sex <- 
  factor(data_tab1$sex, levels=c("1","2", "0"),
         labels=c("Male", 
                  "Female", "Unknown"))

data_tab1$dual <- 
  factor(data_tab1$dual, levels=c("1", "0"),
         labels=c("Eligible", 
                  "Ineligible"))

data_tab1$agerange = as.factor(ifelse(data_tab1$age_dob %in% c(65:69),1, 
                                  ifelse(data_tab1$age_dob %in% c(70:74),2, 
                                         ifelse(data_tab1$age_dob %in% c(75:79),3,  
                                                ifelse(data_tab1$age_dob %in% c(80:84),4, 
                                                       ifelse(data_tab1$age_dob %in% c(85:89),5,
                                                              ifelse(data_tab1$age_dob %in% c(90:94),6,7)))))))

data_tab1$agerange <- 
  factor(data_tab1$agerange, levels=c("1","2","3","4","5","6","7"),
         labels=c("65 - 69", 
                  "70 - 74", 
                  "75 - 79", 
                  "80 - 84", 
                  "85 - 89", 
                  "90 - 94", 
                  "95 + " ))

data_tab1$race <- 
  factor(data_tab1$race_rti, levels=c("1","2","3","4","5","6","0"),
         labels=c("White", 
                  "Black", 
                  "Other", 
                  "Asian", 
                  "Hispanic", 
                  "North American Native", 
                  "Unknown" ))

label(data_tab1$sex)       <- "Sex"
label(data_tab1$race)       <- "Race"
label(data_tab1$agerange)       <- "Age"



units(data_tab1$age_dob)       <- "years"
units(data_tab1$agerange)       <- "years"


caption  <- "Table 1"



my.render.cont <- function(x) {
  with(stats.apply.rounding(stats.default(x), digits=2), c("",
                                                           "Mean (SD)"=sprintf("%s (&plusmn; %s)", MEAN, SD)))
}
my.render.cat <- function(x) {
  c("", sapply(stats.default(x), function(y) with(y,
                                                  sprintf("%d (%0.1f %%)", FREQ, PCT))))
}

my.table1<-table1(~ sex + agerange + race + dual +age_dob,
                  data=data_tab1,
                  overall=c(left="Total"),
                  caption=caption,render.continuous=my.render.cont,
                  render.categorical=my.render.cat
)
