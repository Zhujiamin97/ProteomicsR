#'@title  hela_QC_maxquant
#'@author Jiamin Zhu
#'@description
#' file.dir:txt file path
#'


hela_QC_maxquant <- function(){
  
  file.path <- file.choose()
  
  file.path <- paste0(strsplit(file.path,"txt")[[1]][1],"txt")
  
  files <- list.files(files.dir,full.names = TRUE)
  
  ######################################
  idx1 <- grep("proteinGroups",files)
  proteinGroups <- read.table(files[idx1],header = TRUE,sep = "\t",stringsAsFactors = FALSE)
  mv_idx <- which(proteinGroups$Reverse == "+" | proteinGroups$Potential.contaminant == "+")
  if(length(mv_idx)>0){
    mv_con_res <- proteinGroups[-mv_idx,]
  }else{
    mv_con_res <- proteinGroups
  }
  protein_num <- nrow(mv_con_res)
  
  ######################################
  idx2 <- grep("summary",files)
  summary <- read.table(files[idx2],header = TRUE,sep = "\t",stringsAsFactors = FALSE)
  pep_num <- summary$Peptide.Sequences.Identified[1]
  PSMs <- summary$MS.MS.Identified[1]
  MS2 <- summary$Peaks.Sequenced[1]
  
  ######################################
  idx3 <- grep("evidence",files)
  evidence <- read.table(files[idx3],header = TRUE,sep = "\t",stringsAsFactors = FALSE) 
  mass_error_ppm <- mean(evidence$Uncalibrated...Calibrated.m.z..ppm.)
  
  df <- data.frame(Proteins = protein_num,
                   Peptides = pep_num,
                   PSMs = PSMs,
                   "MS/MS" = MS2,
                   "Mass error (ppm)" = mass_error_ppm)
  return(df)
}