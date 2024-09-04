mult_DIANN_DA_depth <- function(){
  
  pg.tsv.paths <- dir(pattern = "report.pg_matrix.tsv",recursive = TRUE,full.names = TRUE)
  
  pr.tsv.paths <- dir(pattern = "report.pr_matrix.tsv",recursive = TRUE,full.names = TRUE)
  
  Results <- 
    
    lapply(1:length(pg.tsv.paths), function(x){
      
      pg.tsv <- data.table::fread(pg.tsv.paths[x])
      pg.path <- dirname(pg.tsv.paths[x])
      
      pr.tsv <- data.table::fread(pr.tsv.paths[x])
      pr.path <- dirname(pr.tsv.paths[x])
      
      if(pg.path != pr.path){
        stop("项目名不一致")
      }
      
      # report.pg_matrix
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
      Stripped.Sequence <- pr.tsv$Stripped.Sequence %>% unique() %>% length()
      message("---------------------------------------------------")
      print(protein)
      
      print(samples.num)
      
      da <- data.frame(info = "Stripped.Sequence",
                       value = Stripped.Sequence)
      print(da)
      message("---------------------------------------------------")
      
      
      t <- data.frame(info = "--------------------------------------",
                      value = "-----------------")
      
      rbind(protein,samples.num,da,t)
      
    }) %>% do.call(rbind,.)
  
  openxlsx::write.xlsx(Results,"depth.xlsx") 
  message("done!")
}
