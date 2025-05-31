# The Geography of Memory: Regional and Demographic Variabilities in Dementia-Related Diagnostic Code Utilization in U.S. Medicare Hospitalizations
Code for _Federica Spoto, Jiazi Tian, Jonas Hüegel, Daniel Tran Ortega, Christine S. Ritchie, Deborah Blacker, Francesca Dominici, Chirag Patel, Daniel Mork, Hossein Estiri_ "The Geography of Memory: Regional and Demographic Variabilities in Dementia-Related Diagnostic Code Utilization in U.S. Medicare Hospitalizations"

## Abstract


## Workflow 
- analysis: contains the R code for the analyses.
    - 01_data_building_tab1.R: creates the cohort of interest and Table 1.
    - 02_descriptive.R: creates Figure 1 and Figure 2.
    - 03_national_ethnicity_hup.R: computes and compares National HUP with stratification by ethnicity.
    - 04_state_similarity.R: computes HUP by state and compares to National HUP. Plots the result.
    - 05_county_similarity.R: computes HUP by county and compares to National HUP. Plots the result.
    - 06_county_regression.R: prepares the variables and carries out the fixed effect regression analysis.
    
- lib: contains auxiliary functions.
    - tSPMPlus_function.R: runs the tSPM+ algorithm.
    - get_data_graph.R: gets the datagraph.
    - create_corr_matrix.R: creates the correlation matrix from the output of the tSPM+ algorithm.
    - get_dendrogram.R: gets the dendrogram.

