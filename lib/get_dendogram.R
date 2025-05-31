if(!require(pacman)) install.packages("pacman")
require(dplyr)
require(tidyr)
#devtools::install_github("clai-group/mlho")
require(mlho)
require(DT)
require(igraph)
require(ggraph)
require(RColorBrewer)
pacman::p_load(data.table, devtools, backports, Hmisc, tidyr,dplyr,ggplot2,plyr,scales,readr,RcppParallel,
               httr, DT, lubridate, tidyverse,reshape2,foreach,doParallel,caret,gbm,lubridate,praznik,epitools,tcltk)

get_dendogram = function(data_graph){
  
  connect <- data_graph %>% dplyr::select (start, end, mean_val, freq)
  colnames(connect) <- c('from','to','value', 'freq') 
  
  codes <- unique(c(unique(data_graph$start), unique(data_graph$end)))
  
  phen_label <- df %>% dplyr::select (phenotype, CONCEPT_CD) %>% 
    distinct() %>% 
    group_by(CONCEPT_CD) %>%
    dplyr::mutate(phenotype = ifelse(n() > 1, paste(phenotype, collapse = "&"), phenotype)) %>%
    distinct()
  unique(phen_label$phenotype) # 4
  
  phen_label <- subset(phen_label, phen_label$CONCEPT_CD %in% codes)
  
  d1 <- data.frame(from="origin", to=unique(phen_label$phenotype))
  d2 <- phen_label
  colnames(d2)<- c("from", "to")
  hierarchy <- rbind(d1, d2)
  
  
  c(as.character(connect$from), as.character(connect$to)) %>%
    as.tibble() %>%
    group_by(value) %>%
    dplyr::summarize(n=n()) -> nodes
  
  phe <- phen_label  %>% 
    group_by(phenotype) %>%  dplyr::summarize(n=n())
  
  colnames(nodes) <- c("name", "value")
  colnames(phe) <- c("name", "value")
  orign <- data.frame(name ="origin", value="68")
  
  vertices <- rbind(phe, nodes, orign)
  vertices$group  <-  hierarchy$from[match( vertices$name, hierarchy$to ) ]
  vertices$value = as.numeric(vertices$value)
  vertices <- data.frame(vertices)
  
  mygraph <- graph_from_data_frame(hierarchy, vertices=vertices)
  
  from  <-  match(connect$from, vertices$name)
  to  <-  match(connect$to, vertices$name)
  
  
  p <- ggraph(mygraph, layout = 'dendrogram', circular = TRUE) +
    geom_conn_bundle(data = get_con(from = from, to = to), aes(colour= connect$value[after_stat(index)])) +
    scale_edge_color_continuous(low="white", high="red")
  #scale_edge_colour_distiller(palette = "RdPu")
  # width = connect$freq[after_stat(index)]
  
  
  p = p+ 
    geom_node_point(aes(filter = leaf, x = x*1.05, y=y*1.05, colour=group, size=vertices$value)) +
    scale_colour_manual(values= rep( brewer.pal(9,"Paired") , 30)) +
    scale_size_continuous(range = c(0.1,10) )  +
    geom_node_text(aes(label = name),  colour = 'black', size=3, nudge_x = p$data$x * .15, nudge_y = p$data$y * .15) +
    theme_void()
  
  
  return(p)
}
