#'@title Specrtonaut_DA_depth
#'@data 2024/09/03
#'

Spectronaut_DA_depth <- function(){
  
  tsv.path <- file.choose()
  
  # workpath <- dirname(tsv.path)
  
  quant_protein <- data.table::fread(tsv.path)
  
  PG.ProteinGroups <- quant_protein$PG.ProteinGroups %>% length()
  
  samples <- quant_protein[,-1] %>% as.data.frame()
  samples.num <- lapply(1:ncol(samples),function(x){
    
    info <- colnames(samples)[x]
    info <- strsplit(info,"] ",fixed = TRUE)[[1]] %>% tail(.,1)
    value <- sum(!is.na(samples[,x])) 
    data.frame(info,value)
    
  }) %>% do.call(rbind,.)
 
  message("----------------------------------------------------------------------")
  message("PG.ProteinGroups:")
  print(PG.ProteinGroups)
  message("----------------------------------------------------------------------")
  print(samples.num)
  message("----------------------------------------------------------------------")
  message("done!")
}


#'@title DA_depth
#'@data 2024/09/03
#'
DA_depth <- function(){
  
  tsv.path <- file.choose()
  
  if(tsv.path %like% "Protein_Quant.tsv"){
    message("----------------------------------------------------------------------")
    message("Spectronaut Results:")
    
    quant_protein <- data.table::fread(tsv.path)
    
    PG.ProteinGroups <- quant_protein$PG.ProteinGroups %>% length()
    
    samples <- quant_protein[,-1] %>% as.data.frame()
    samples.num <- lapply(1:ncol(samples),function(x){
      
      info <- colnames(samples)[x]
      info <- strsplit(info,"] ",fixed = TRUE)[[1]] %>% tail(.,1)
      value <- sum(!is.na(samples[,x])) 
      data.frame(info,value)
      
    }) %>% do.call(rbind,.)
    
    message("----------------------------------------------------------------------")
    message("PG.ProteinGroups:")
    print(PG.ProteinGroups)
    message("----------------------------------------------------------------------")
    print(samples.num)
    message("----------------------------------------------------------------------")
    message("done!")
  }else{
    message("---------------------------------------------------")
    message("DIA-NN Results:")
    
    workpath <- dirname(tsv.path)
    
    # report.pg_matrix
    pg.tsv <- data.table::fread(paste0(workpath,"/report.pg_matrix.tsv"))
    protein <- data.frame(info = "protein",
                          value =  pg.tsv$Protein.Group %>% length())
    samples <- pg.tsv[,-c(1:5)] %>% as.data.frame()
    samples.num <- lapply(1:ncol(samples),function(x){
      
      sample.path <- colnames(samples)[x]
      info <- strsplit( sample.path,"\\",fixed = TRUE)[[1]] %>% tail(.,1)
      value <- sum(!is.na(samples[,x])) 
      data.frame(info,value)
      
    }) %>% do.call(rbind,.)
    # report.pg_matrix
    pr.tsv <- data.table::fread(paste0(workpath,"/report.pr_matrix.tsv"))
    Stripped.Sequence <- pr.tsv$Stripped.Sequence %>% unique() %>% length()
    message("---------------------------------------------------")
    print(protein)
    message("---------------------------------------------------")
    print(samples.num)
    message("---------------------------------------------------")
    print(data.frame(info = "Stripped.Sequence",
                     value = Stripped.Sequence))
    message("---------------------------------------------------")
    message("done!")
    
  }
  
}