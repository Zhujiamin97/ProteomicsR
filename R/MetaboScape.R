#'@title metabo

MetaboScape_select_name <- function(){
  
  # 样本信息文件路径
  filepath <- file.choose()
  
  workpath <- dirname(filepath)
  
  setwd(workpath)
  
  list <- openxlsx::read.xlsx(filepath,sheet = "Sheet1",colNames = FALSE) %>% as.data.frame()
  
  #negative
  if(file.exists("neg_quant.xlsx")){
    
    neg_quant <- read_xlsx("neg_quant.xlsx")
    
    if(nrow(list) != ncol(neg_quant[,-c(1:12)])){
      stop("请检查列表名称是否数量相等！！")
    }
    
    neg_quant <- cbind(neg_quant[,1:12],select(neg_quant,list[,1]))
    openxlsx::write.xlsx(neg_quant,"neg_quant_select.xlsx")
    message("-------------------------")
    message(paste0("Negative模式鉴定数量:",length(unique(neg_quant$Name))))
    message("-------------------------")
    
  }else{
    message("no neg_quant.xlsx file!")
  }
  # positive
  if(file.exists("pos_quant.xlsx")){
    
    pos_quant <- read_xlsx("pos_quant.xlsx")
    
    if(nrow(list) != ncol(pos_quant[,-c(1:12)])){
      stop("请检查列表名称是否数量相等！！")
    }
    
    pos_quant <- cbind(pos_quant[,1:12],select(pos_quant,list[,1]))
    openxlsx::write.xlsx(pos_quant,"pos_quant_select.xlsx")
    message(paste0("Positive模式鉴定数量:",length(unique(pos_quant$Name))))
    message("-------------------------")
  }else{
    message("no pos_quant.xlsx file!")
    
  }
  message("DONE!")
  
}