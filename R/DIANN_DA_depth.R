#'@title diann_DA_depth
#'@data 2024/09/02
#'

diann_DA_depth <- function(){
  
  tsv.path <- file.choose()
  
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