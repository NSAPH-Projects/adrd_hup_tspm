library(dplyr)
library(tidyverse)
library(magrittr)
library(lubridate)
require(tidyr)

create_corr_matrix = function(data_graph, all_codes = c("3310",  "29420", "33119", "33183", "33182", "2900", "29021", "29421",
                                                        "797", "2903",  "3314", "29042", "29410", "29043", "2948", "29020", 
                                                        "29010", "29041", "29411", "29040", "3319", "3315", "2949", "2909", 
                                                        "29011", "33189", "3313", "29012", "33111", "3316", "2940", "3312")){
  matrix_checkentry = data_graph[,c("start", "end", "mean_val")] %>%
    mutate(entry = paste0(start,"-",end))
  
  for (start_code in all_codes){
    for (end_code in all_codes){
      if (! paste0(start_code,"-",end_code) %in% matrix_checkentry$entry){
        matrix_checkentry = matrix_checkentry %>% 
          add_row(start = (start_code), end = end_code, mean_val = 0, entry = paste0(start_code,"-",end_code))
      }
    }
  }
  
  corr = acast(matrix_checkentry, start~end, value.var="mean_val", mean)
  return(corr)
}