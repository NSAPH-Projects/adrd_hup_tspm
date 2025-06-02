# Load libraries
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

setwd("/n/dominici_nsaph_l3/Lab/projects/adrd_hup_tspm/analysis")

source("../lib/get_data_graph.R")
source("../lib/get_dendogram.R")
source("../lib/create_corr_matrix.R")

load("../data/scratch/ADRD_hospitalization_2016-2018_clin.RData")

# Prepare data ------------------------------------------------------------

# Import shapefile
shapefile_path <- "/n/dominici_nsaph_l3/Lab/projects/adrd_hup_tspm/data/input/tl_2019_us_county/tl_2019_us_county.shp"
county_data <- st_read(shapefile_path) %>% 
  #filter(!STATEFP %in% c( "02","15", "66", "72", "60", "69", "78")) %>%
  dplyr::select(c("STATEFP", "GEOID", "geometry"))

# Load denom files to get counties
files_denom <- c("/n/dominici_nsaph_l3/Lab/projects/analytic/mbsf_medpar_denom/mbsf_medpar_denom_2016.parquet",
                 "/n/dominici_nsaph_l3/Lab/projects/analytic/mbsf_medpar_denom/mbsf_medpar_denom_2017.parquet",
                 "/n/dominici_nsaph_l3/Lab/projects/analytic/mbsf_medpar_denom/mbsf_medpar_denom_2018.parquet")
denom <- as.data.frame(rbindlist(lapply(files_denom,
                                        read_parquet)))

fips_df = denom %>%
  select(zip, county, year)%>%
  group_by(zip, year)%>% 
  filter(row_number() == 1) %>%
  ungroup()%>%
  na.omit()

# Get hospitalizations' county
ADRD_hospitalization = left_join(ADRD_hospitalization, fips_df)
fips_ADRD = unique(ADRD_hospitalization$county)
fips = unique(county_data$GEOID)
valid_fips = fips_ADRD[fips_ADRD %in% fips]
ADRD_hospitalization = ADRD_hospitalization[ADRD_hospitalization$county %in% valid_fips,]

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


# Compute and compare HUP -------------------------------------------------------------

# National
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

# Create dataframe to store results
county_sk = data.frame(county = valid_fips,
                       random_skewers = NA,
                       mad = NA,
                       n = NA)
# Compute and compare HUP per every county
for(county in valid_fips){
  gc()
  # Filter on county
  ADRD_county = ADRD_hospitalization[ADRD_hospitalization$county == county,]
  
  # Numper of patients per county
  county_sk[ county_sk$county == county , 'n'] = ADRD_county %>%
    summarize(n = n_distinct(bene_id))
  
  # Censor county with less than 11 patients
  if(county_sk[ county_sk$county == county , 'n'] >10){
    dbmart_county = ADRD_county %>%
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
    
    data_graph_county = get_data_graph(dbmart_county, df)
    
    # Create correlation matrix
    corr_county = tryCatch({create_corr_matrix(data_graph_county, df$CONCEPT_CD)}, error =  function(err){f <- NA
    return(f)})
    
    # Compute similarity
    county_sk[ county_sk$county == county , 'random_skewers'] = tryCatch({skewers(corr_county, corr_us, nsim=1000, method=NULL)$r}, error =  function(err){f <- NA; return(f)})
    
    county_sk[ county_sk$county == county , 'mad'] = tryCatch({mean(abs(corr_county - corr_us))}, error =  function(err){f <- NA; return(f)})
  }
  
  save(county_sk, file = "../data/county_corr_to_national_clin_groups.RData")
}

save(county_sk, file = "../data/county_corr_to_national_clin_groups.RData")

# Results -----------------------------------------------------------------
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

county_sk_sh = county_sk_sh %>%
  mutate(mad_compl = 1 - mad)

## Plots -----

# Random Skewers
county_rs = county_sk_sh %>% 
  st_as_sf() %>% 
  st_simplify() %>% 
  ggplot() + 
  geom_sf(aes(fill = random_skewers), linewidth = 0.1, col = "black") +  
  scale_fill_viridis_c( "Similarity ", na.value = "lightgrey", limits = c(0,1))+
  theme_void()+ 
  ggtitle("(b) National HUP vs County HUP") + 
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
ggsave("county_rs_5groups.png", county_rs, width = 10, height = 8)
ggsave("county_rs_5groups.pdf", county_rs, width = 10, height = 8, device = cairo_pdf)

county_sk_sh %>% 
  st_as_sf() %>% 
  st_simplify() %>% 
  ggplot() + 
  geom_sf(aes(fill = mad), linewidth = 0.1, col = "black") +  
  scale_fill_viridis_c( "MAD", na.value = "lightgrey")+
  theme_void()+ 
  ggtitle("Mean absolute difference - National vs County") + 
  theme(plot.title = element_text(size = 22, hjust = 0.5),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        line = element_blank(),
        axis.title = element_blank(),
        legend.text = element_text( size = 16),
        legend.title = element_text(size = 18),
        legend.key.width = unit(50, "points"),
        panel.grid.major = element_line(colour = "transparent"), 
        legend.position = "bottom",
        legend.direction = "horizontal") 

county_mad = county_sk_sh %>% 
  st_as_sf() %>% 
  st_simplify() %>% 
  ggplot() + 
  geom_sf(aes(fill = mad_compl), linewidth = 0.1, col = "black") +  
  scale_fill_viridis_c( "Similarity ", na.value = "lightgrey", breaks = c(0.75, 0.85, 0.95))+
  theme_void()+ 
  ggtitle("(b) National HUP vs County HUP") + 
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
ggsave("county_mad_5groups.png", county_mad, width = 10, height = 8) 
ggsave("county_mad_5groups.pdf", county_mad, width = 10, height = 8, device = cairo_pdf) 


# Association between MAD and RS
ggplot(county_sk_sh, aes(x = mad_compl, y = random_skewers)) +
  geom_point() +
  labs(
    title = "Scatter Plot of mad_compl vs. random_skewers",
    x = "MAD Completeness",
    y = "Random Skewers"
  ) +
  theme_minimal()


# Numerosity
summary(county_sk_sh$n)

ggplot(county_sk_sh, aes(x = n)) +
  geom_histogram(binwidth = 100, fill = "skyblue", color = "black") +
  labs(title = "Patients per county", x = "Number of patients", y = "Frequency") +
  theme_minimal()


