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

setwd("/n/dominici_nsaph_l3/Lab/projects/adrd_hup_tspm/analysis")

source("../lib/get_data_graph.R")
source("../lib/get_dendogram.R")
source("../lib/create_corr_matrix.R")

load("../data/scratch/ADRD_hospitalization_2016-2018_clin.RData")

# Get ADRD codes 
dementia_codes <- read.csv("../data/ADRD_ICD10.csv", stringsAsFactors = FALSE)
names(dementia_codes) <- c("ICD10", "category", "Description")
dementia_codes$ICD10[1] = "F0150"
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
corr_us = create_corr_matrix(data_graph_us, df$CONCEPT_CD)
save(corr_us, file = "../data/corr_us_clin_groups.RData")

# States HUP 
states_sf <- states( cb=TRUE) %>% filter(!STATEFP %in% c( "02","15", "66", "72", "60", "69", "78"))
states_sf[ , 'state_us_skewers'] = NA
states_sf[ , 'state_us_mad'] = NA

states_sf %<>%
  dplyr::select(c("STUSPS", "NAME", "state_us_skewers", "state_us_mad", "geometry"))

states = unique(states_sf$STUSPS)

for(state in states){
  ADRD_state = ADRD_hospitalization[ADRD_hospitalization$state == state,]
  
  # State ----
  dbmart_state = ADRD_state %>%
    ungroup()%>%
    dplyr::select("bene_id", "admission_date", paste0("DIAG",1:25))%>%
    pivot_longer(cols = starts_with('DIAG'))%>%
    filter(value %in% ADRD_codes)%>%
    dplyr::rename(patient_num = bene_id,
                  start_date = admission_date,
                  phenx = value)%>%
    left_join(dementia_groups) %>%
    dplyr::select(-phenx) %>%
    dplyr::rename(phenx = group)
  
  data_graph_state = get_data_graph(dbmart_state, df)
  
  corr_state = create_corr_matrix(data_graph_state, df$CONCEPT_CD)
  
  # Similarity to national HUP
  states_sf[states_sf$STUSPS == state , 'state_us_skewers'] = skewers(corr_us, corr_state, nsim=1000, method=NULL)$r
  states_sf[states_sf$STUSPS == state , 'state_us_mad'] = mean(abs(corr_state - corr_us))
  
  save(states_sf, file = "../data/states_corr_clin_groups.RData")
}

load( "../data/states_corr_clin_groups.RData")

# Complement of MAD
states_sf = states_sf%>%
  mutate(mad_compl = 1 - state_us_mad)%>%
  st_transform(crs = 5070)

# Plot results
state_rs = states_sf %>% 
  st_as_sf() %>% 
  st_simplify() %>% 
  ggplot() + 
  geom_sf(aes(fill = state_us_skewers), linewidth = 0.1, col = "black") +  
  scale_fill_viridis_c( "Similarity ", na.value = "lightgrey",  limits = c(0, 1) )+
  theme_void()+ 
  ggtitle("(a) National HUP vs State HUP") + 
  theme(plot.title = element_text(size = 30, hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        line = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text( size = 22),
        legend.title = element_text(size = 25),
        legend.key.width = unit(50, "points"),
        panel.grid.major = element_line(colour = "transparent"), 
        legend.position = "bottom",
        legend.direction = "horizontal") 

ggsave("state_rs_5groups.png", state_rs, width = 10, height = 8)
ggsave("state_rs_5groups.pdf", state_rs, width = 10, height = 8, device = cairo_pdf)

state_mad = states_sf %>% 
  st_as_sf() %>% 
  st_simplify() %>% 
  ggplot() + 
  geom_sf(aes(fill = mad_compl), linewidth = 0.1, col = "black") +  
  scale_fill_viridis_c( "Similarity ", na.value = "lightgrey", breaks = c(0.94, 0.96, 0.98),labels = c("0.94", "0.96", "0.98"))+
  theme_void()+ 
  ggtitle("(a) National HUP vs State HUP") + 
  theme(plot.title = element_text(size = 30, hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        line = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text( size = 22),
        legend.title = element_text(size = 25),
        legend.key.width = unit(50, "points"),
        panel.grid.major = element_line(colour = "transparent"), 
        legend.position = "bottom",
        legend.direction = "horizontal") 

ggsave("state_mad_5groups.png", state_mad, width = 10, height = 8)
ggsave("state_mad_5groups.pdf", state_mad, width = 10, height = 8, device = cairo_pdf)


