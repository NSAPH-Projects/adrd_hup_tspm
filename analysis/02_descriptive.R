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
library(tigris)
library(sf)
library(corrplot)

setwd("/n/dominici_nsaph_l3/Lab/projects/adrd_hup_tspm/analysis")

source("../lib/get_data_graph.R")
source("../lib/get_dendogram.R")
source("../lib/create_corr_matrix.R")
source("../lib/tSPMPlus_function.R")

load("../data/scratch/ADRD_hospitalization_2016-2018_clin.RData")
states_sf <- states( cb=TRUE) %>% filter(!STATEFP %in% c( "02","15", "66", "72", "60", "69", "78")) %>%
  dplyr::select(c("STUSPS", "NAME", "geometry"))
states = unique(states_sf$STUSPS)

# Read ADRD code lookup
dementia_codes <- read.csv("../data/ADRD_ICD10.csv", stringsAsFactors = FALSE)
names(dementia_codes) <- c("ICD10", "category", "Description")
ADRD_codes = dementia_codes$ICD10

# Figure 1: Frequency of ICD-10 codes -------------------------------------

# Create a vector with all observed ICD-10 codes in hospitalizations
diags = as.vector(as.matrix(ADRD_hospitalization[,paste0("DIAG", 1:25)]))
diags_ADRD = diags[diags %in% ADRD_codes] #filter on ADRD-related

# Frequences
diag_counts <- as.data.frame(table(diags_ADRD), stringsAsFactors=FALSE) %>%
  dplyr::rename(ICD10 = diags_ADRD,
         Freq  = Freq)

# Merge diagnostic category
diag_counts <- diag_counts %>%
  left_join(dementia_codes %>% 
              select(ICD10, category),
            by = "ICD10")
# Order
diag_counts <- diag_counts %>%
  arrange(category, desc(Freq)) %>%
  mutate(ICD10 = factor(ICD10, levels = unique(ICD10)))

# Plot
ggplot(diag_counts, aes(x = ICD10, y = Freq, fill = category)) +
  geom_col(color = "black") +
  geom_text(aes(label = Freq),
            vjust = -0.3,
            size = 3,
            show.legend = FALSE) +
  scale_fill_brewer(palette = "Set2", na.value = "gray80") +
  theme_minimal() +
  theme(
    axis.text.x      = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10),
    axis.text.y      = element_text(size = 10),
    axis.title       = element_text(face = "bold", size = 12),
    plot.title       = element_text(face = "bold", size = 14, hjust = 0.5)
  ) +
  labs(
    x     = "ICD-10 code",
    y     = "# unique patients with claim",
    fill  = "Diagnostic category",
    title = "ICD-10 codes frequency in hospitalization claims"
  )

# Figure 2. State percentage of diagnosis category compared to the nation --------
# Create df to input in tSPM+
dementia_groups = dementia_codes %>%
  select(ICD10, category) %>%
  dplyr::rename(phenx = ICD10) %>%
  mutate(group = category)

df <- dementia_groups %>%
  dplyr::select(c("category", "group")) %>% 
  dplyr::rename(phenotype = group, CONCEPT_CD = category)
df = unique(df)


# National ----
dbmart = ADRD_hospitalization %>%
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

### TSPM+ -----------------
tspm_res = tSPMPlus_function(dbmart, df)
corseq = tspm_res$corseq
corrs_final= tspm_res$corrs_final
phenxlookup = tspm_res$phenxlookup %>%
  mutate(num_Phenx = as.character(num_Phenx))

# Create a data frame connect that includes the start and end phenX of each sequence
corseq <- as.data.frame(corseq)
corseq$patient_num <- as.character(corseq$patient_num)
corseq$sequence <- as.character(corseq$sequence)
num_pt <- n_distinct(corseq$patient_num)
nodes <- corseq %>% dplyr::select( startPhen, endPhenx, patient_num) %>%
  distinct() %>%
  gather(key = "node_type", value = "num_Phenx", startPhen, endPhenx) %>%
  group_by(num_Phenx) %>%
  dplyr::summarize(total_frequency = n_distinct(patient_num)) %>%
  mutate(national = total_frequency/num_pt* 100) %>%
  left_join(phenxlookup) %>%
  select(c(phenx, national))

# Create dataframe to store results
res = phenxlookup %>%
  select(phenx) %>%
  left_join(nodes)

for(state in states){
  
  # State ----
  ADRD_state = ADRD_hospitalization[ADRD_hospitalization$state == state,]
  
  dbmart_state = ADRD_state %>%
    ungroup()%>%
    dplyr::select("bene_id", "admission_date", paste0("DIAG",1:25))%>%
    pivot_longer(cols = starts_with('DIAG'))%>%
    filter(value %in% ADRD_codes)%>%
    dplyr::rename(patient_num = bene_id,
                  start_date = admission_date,
                  phenx = value) %>%
    left_join(dementia_groups) %>%
    dplyr::select(-phenx) %>%
    dplyr::rename(phenx = clin)
  tspm_res_state = tSPMPlus_function(dbmart_state, df)
  corseq_state = tspm_res_state$corseq
  corrs_final_state= tspm_res_state$corrs_final
  phenxlookup_state = tspm_res_state$phenxlookup%>%
    mutate(num_Phenx = as.character(num_Phenx))
  
  # Create a data frame connect that includes the start and end phenX of each sequence
  corseq_state <- as.data.frame(corseq_state)
  corseq_state$patient_num <- as.character(corseq_state$patient_num)
  corseq_state$sequence <- as.character(corseq_state$sequence)
  
  num_pt_state <- n_distinct(corseq_state$patient_num)
  
  nodes <- corseq_state %>% dplyr::select( startPhen, endPhenx, patient_num) %>%
    distinct() %>%
    gather(key = "node_type", value = "num_Phenx", startPhen, endPhenx) %>%
    group_by(num_Phenx) %>%
    dplyr::summarize(total_frequency = n_distinct(patient_num)) %>%
    mutate(state = total_frequency/num_pt_state* 100) %>%
    left_join(phenxlookup_state) %>%
    select(c(phenx, state)) %>%
    rename_with(~c(state), c(state))
  
  # store results
  res = res %>%
    left_join(nodes)
}

res[is.na(res)] = 0

data_long <- res %>%
  select(-national) %>%
  pivot_longer(cols = -phenx, names_to = "State", values_to = "Value")

# Plot
p1 = ggplot(data_long, aes(x = phenx, y = Value, color = State)) +
  geom_point(size = 3) +
  geom_segment(
    # note: linewidth, not size
    data    = res,
    aes(
      # convert factor to numeric so we know its x‐position
      x    = as.numeric(phenx) - 0.1,
      xend = as.numeric(phenx) + 0.1,
      y    = national,
      yend = national
    ),
    linewidth = 1.5,   # thickness of the bar
    color     = "black"
  ) +
  scale_x_discrete() +
  scale_color_viridis_d() +
  labs(
    title = "Percentage of patients with at least one hospitalization VS groups",
    x     = "Diagnosis category",
    y     = "Percentage"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("hosp.png", p1, width = 10, height = 6)

data_long %>% dplyr::group_by(as.character(phenx)) %>% 
  dplyr::summarise(min = min(Value), max = max(Value))


