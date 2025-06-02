library(phytools)
library(fst)
library(data.table)
library(dplyr)
library(tidyverse)
library(magrittr)
library(lubridate)
require(tidyr)
library(readr)
library(arrow)
library(corrplot)
library(tigris)
library(ggplot2)
library(sf)
sf_use_s2(FALSE)
# library(raster)

setwd("/n/dominici_nsaph_l3/Lab/projects/adrd_hup_tspm/analysis")

source("../lib/get_data_graph.R")
source("../lib/get_dendogram.R")
source("../lib/create_corr_matrix.R")

load("../data/scratch/ADRD_hospitalization_2016-2018_clin.RData")

# Get ADRD codes 
dementia_codes <- read.csv("../data/ADRD_ICD10.csv", stringsAsFactors = FALSE)
names(dementia_codes) <- c("ICD10", "category", "Description")
ADRD_codes = dementia_codes$ICD10

dementia_groups = dementia_codes %>%
  select(ICD10, category) %>%
  dplyr::rename(phenx = ICD10) %>%
  mutate(group = category)

df <- dementia_groups %>%
  dplyr::select(c("category", "group")) %>% 
  dplyr::rename(phenotype = group, CONCEPT_CD = category)
df = unique(df)

# National HUP 
dbmart_us = ADRD_hospitalization %>%
  ungroup()%>%
  dplyr::select("bene_id", "admission_date", paste0("DIAG",1:25))%>%
  pivot_longer(cols = starts_with('DIAG'))%>%
  filter(value %in% ADRD_codes)%>%
  dplyr::rename(patient_num = bene_id,
                start_date = admission_date,
                phenx = value) %>%
  left_join(dementia_groups) %>%
  dplyr::select(-phenx) %>%
  dplyr::rename(phenx = group)

data_graph_us = get_data_graph(dbmart_us, df)

## Get HUP per ethnicity ------
# Black
ADRD_hospitalization_black = ADRD_hospitalization %>%
  filter(race_rti==2)

dbmart_black = ADRD_hospitalization_black %>%
  ungroup()%>%
  select("bene_id", "admission_date", paste0("DIAG",1:25))%>%
  pivot_longer(cols = starts_with('DIAG'))%>%
  filter(value %in% ADRD_codes)%>%
  dplyr::rename(patient_num = bene_id,
                start_date = admission_date,
                phenx = value) %>%
  left_join(dementia_groups) %>%
  dplyr::select(-phenx) %>%
  dplyr::rename(phenx = group)

data_graph_black = get_data_graph(dbmart_black, df)

# Hispanic
ADRD_hospitalization_hisp = ADRD_hospitalization %>%
  filter(race_rti==5)

dbmart_hisp = ADRD_hospitalization_hisp%>%
  ungroup()%>%
  select("bene_id", "admission_date", paste0("DIAG",1:25))%>%
  pivot_longer(cols = starts_with('DIAG'))%>%
  filter(value %in% ADRD_codes)%>%
  dplyr::rename(patient_num = bene_id,
                start_date = admission_date,
                phenx = value) %>%
  left_join(dementia_groups) %>%
  dplyr::select(-phenx) %>%
  dplyr::rename(phenx = group)

data_graph_hisp = get_data_graph(dbmart_hisp, df)

# White
ADRD_hospitalization_white = ADRD_hospitalization %>%
  filter(race_rti==1)

dbmart_white = ADRD_hospitalization_white%>%
  ungroup()%>%
  select("bene_id", "admission_date", paste0("DIAG",1:25))%>%
  pivot_longer(cols = starts_with('DIAG'))%>%
  filter(value %in% ADRD_codes)%>%
  dplyr::rename(patient_num = bene_id,
                start_date = admission_date,
                phenx = value) %>%
  left_join(dementia_groups) %>%
  dplyr::select(-phenx) %>%
  dplyr::rename(phenx = group)
data_graph_white = get_data_graph(dbmart_white, df)


# Asian
ADRD_hospitalization_asian = ADRD_hospitalization %>%
  filter(race_rti==4)

dbmart_asian = ADRD_hospitalization_asian%>%
  ungroup()%>%
  select("bene_id", "admission_date", paste0("DIAG",1:25))%>%
  pivot_longer(cols = starts_with('DIAG'))%>%
  filter(value %in% ADRD_codes)%>%
  dplyr::rename(patient_num = bene_id,
                start_date = admission_date,
                phenx = value) %>%
  left_join(dementia_groups) %>%
  dplyr::select(-phenx) %>%
  dplyr::rename(phenx = group)

data_graph_asian = get_data_graph(dbmart_asian, df)

## Create correlation matrices -----
corr_us = create_corr_matrix(data_graph_us, df$CONCEPT_CD)
corr_black = create_corr_matrix(data_graph_black, df$CONCEPT_CD)
corr_hisp = create_corr_matrix(data_graph_hisp, df$CONCEPT_CD)
corr_white = create_corr_matrix(data_graph_white, df$CONCEPT_CD)
corr_asian = create_corr_matrix(data_graph_asian, df$CONCEPT_CD)

which.max(corr_us)
data_graph_us[order(data_graph_us$mean_val,decreasing = TRUE),]

corrplot(corr_us, title = "National Correlations", mar=c(0,0,2,0), method = 'color', addCoef.col = "black", number.digits = 2)
corrplot(corr_white, title = "White cohort", mar=c(0,0,2,0), method = 'color', addCoef.col = "black", number.digits = 2)
corrplot(corr_black, title = "Black cohort", mar=c(0,0,2,0), method = 'color', addCoef.col = "black", number.digits = 2)
corrplot(corr_hisp, title = "Hispanic cohort", mar=c(0,0,2,0), method = 'color', addCoef.col = "black", number.digits = 2)
corrplot(corr_asian, title = "Asian cohort", mar=c(0,0,2,0), method = 'color', addCoef.col = "black", number.digits = 2)

# Number of hospitalizations
length(unique(ADRD_hospitalization_white$bene_id))
nrow(ADRD_hospitalization_white)

length(unique(ADRD_hospitalization_black$bene_id))
nrow(ADRD_hospitalization_black)

length(unique(ADRD_hospitalization_hisp$bene_id))
nrow(ADRD_hospitalization_hisp)

length(unique(ADRD_hospitalization_asian$bene_id))
nrow(ADRD_hospitalization_asian)

# Similarity with national HUP
skewers(corr_white, corr_us, nsim=1000, method=NULL)$r
skewers(corr_black, corr_us, nsim=1000, method=NULL)$r
skewers(corr_hisp, corr_us, nsim=1000, method=NULL)$r
skewers(corr_asian, corr_us, nsim=1000, method=NULL)$r
