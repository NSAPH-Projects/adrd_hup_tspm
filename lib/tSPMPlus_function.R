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

tSPMPlus_function = function(dbmart, df){
  
  db <- tSPMPlus::transformDbMartToNumeric(dbmart)
  phenxlookup <- db$phenxLookUp
  patlookup <- db$patientLookUp
  
  phenxOfInterest = c(as.numeric(db$phenxLookUp$num_Phenx)) ### for the dementia/AD work we only are interested in those codes. 
  ## it should be other codes of interest and not general at this point
  temporalBucket =  c(0,1,3)
  minDuration = 0 #techical parameter, ignore for now ##TODO if not working with 0 use 1
  bitShift = 0 #techical parameter, ignore for now
  lengthOfPhenx = 7 #techical parameter
  storeSequencesDuringCreation = FALSE  #if true, old way -> writing out "plain" sparse sequences in patient based files, FALSE-> do in memory sparsity
  
  mem_buffer <- 1 #in GB. just a buffer to make sure the computer wont crash
  cores_buffer <- 1
  buffer<- 100000000*mem_buffer
  dbmart_adapt <- tSPMPlus::splitdbMartInChunks(db$dbMart, includeCorSeq = TRUE, buffer = buffer)
  numOfChunks = length(dbmart_adapt$chunks)
  dbmart_num <- dbmart_adapt$chunks[[1]]
  
  sparsity = 0.0001
  numOfThreads = detectCores()-cores_buffer
  
  corseq <- tSPMPlus::getCandidateSequencesForPOI(dbmart_num,
                                                  minDuration,
                                                  bitShift,
                                                  lengthOfPhenx,
                                                  temporalBucket,
                                                  phenxOfInterest,
                                                  storeSequencesDuringCreation,
                                                  numOfThreads = numOfThreads,
                                                  sparsityValue = sparsity)
  
  
  corseq <- dplyr::distinct(corseq, .keep_all = TRUE)
  
  
  corseq <- corseq %>%
    dplyr::group_by(patient_num,sequence,endPhenx,durationBucket) %>%
    dplyr::summarise(count=length((patient_num)))
  
  corseq$value.var <- 1
  
  corseq$startPhen <-substr(corseq$sequence, start = 1, stop = (nchar(corseq$sequence)-lengthOfPhenx))
  
  corseq$startPhen_dur <- paste0(corseq$startPhen,"-",corseq$durationBucket)
  
  dat <- data.table(corseq[c("patient_num","endPhenx","value.var","startPhen_dur")])
  wide.dat.start <- dcast.data.table(dat, patient_num ~ startPhen_dur, value.var="value.var", fun=length)
  
  end <- phenxOfInterest
  
  corseq$value.var <- 1
  
  corseq$startPhen <-substr(corseq$sequence, start = 1, stop = (nchar(corseq$sequence)-lengthOfPhenx))
  
  corseq$startPhen_dur <- paste0(corseq$startPhen,"-",corseq$durationBucket)
  
  dat <- data.table(corseq[c("patient_num","endPhenx","value.var","startPhen_dur")])
  wide.dat.start <- dcast.data.table(dat, patient_num ~ startPhen_dur, value.var="value.var", fun=length)
  
  end <- phenxOfInterest
  
  for (i in seq(1:numOfChunks)) {
    gc()
    cores<-detectCores()
    cl <- parallel::makeCluster(cores) #parallel::makeCluster(cores-cores_buffer)
    doParallel::registerDoParallel(cl)
    tryCatch({
      corrs <- foreach(j = 1:length(end),#
                       .combine='rbind',
                       .multicombine=TRUE,
                       .packages = c("data.table")) %dopar% {
                         tryCatch({
                           dat.j <- data.table(corseq[corseq$endPhenx %in% end[j],c("patient_num","endPhenx","value.var")])
                           lab.j <- dcast.data.table(dat.j, patient_num ~ endPhenx, value.var="value.var", fun=length)
                           colnames(lab.j)[2] <- "label"
                           wide.j <- merge(wide.dat.start,lab.j,by="patient_num",all.x = T)
                           wide.j[is.na(wide.j)] <- 0
                           wide.j$patient_num <- NULL
                           
                           
                           out <- apply(wide.j[, -("label")], 2, cor.test, wide.j$label, method="spearman")
                           p <- data.frame(sapply(out, "[[", "p.value"))
                           p$startPhen_dur <- rownames(p)
                           rownames(p) <- NULL
                           colnames(p)[1] <- "p.value"
                           p$p.adjust <- p.adjust(p$p.value, method = "holm", n = nrow(p))# consider "holm" or "hochberg" or "bonferroni"
                           
                           rho <- data.frame(sapply(out, "[[", "estimate"))
                           rho$startPhen_dur <- rownames(rho)
                           rownames(rho) <- NULL
                           colnames(rho)[1] <- "rho"
                           rho$startPhen_dur <- substr(rho$startPhen_dur,1,nchar(rho$startPhen_dur)-4)
                           cor.j <- merge(rho,p,by="startPhen_dur")
                           cor.j$rho.abs <- abs(cor.j$rho)
                           cor.j$endPhenx <- end[j]
                           
                           
                           rm(rho,p);gc()
                           
                           
                           cor.j
                           
                         },
                         error = function(foll) {cat("ERROR :",conditionMessage(foll), "\n")})
                       }
      
      
      corrs <- merge(corrs,db$phenxLookUp,by.x = "endPhenx",by.y ="num_Phenx" )
      
      corrs$startPhen <- sub("\\-.*", "", corrs$startPhen_dur)
      corrs$sequence <- ifelse (corrs$endPhenx <10,paste0(corrs$startPhen,"000000",corrs$endPhenx),NA)
      corrs$sequence <- ifelse(corrs$endPhenx <100 & corrs$endPhenx >9 ,paste0(corrs$startPhen,"00000",corrs$endPhenx),corrs$sequence)
      corrs$sequence <- ifelse(corrs$endPhenx <1000 & corrs$endPhenx >99 ,paste0(corrs$startPhen,"0000",corrs$endPhenx),corrs$sequence)
      
      corrs$sequence <- as.numeric(corrs$sequence)
      
    },
    error = function(foll) {cat("ERROR in chunk ",i, ": ", conditionMessage(foll), "\n")})
    stopCluster(cl)
    
  }
  
  corrs_final <- corrs %>%
    dplyr::select(-phenx) %>%
    dplyr::group_by(endPhenx,startPhen_dur,sequence,startPhen) %>%
    dplyr::summarise_all(mean, na.rm=TRUE)
  
  return(list(corseq = corseq, corrs_final = corrs_final, phenxlookup = phenxlookup))
}