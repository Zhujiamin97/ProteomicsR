#'@title  hela_QC_maxquant
#'@author Jiamin Zhu
#'@description
#' file.dir:txt file path
#'


hela_QC_maxquant <- function(){
  
  file.path <- file.choose()
  
  file.path <- paste0(strsplit(file.path,"txt")[[1]][1],"txt")
  
  files <- list.files(file.path,full.names = TRUE)
  
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

#'@title Spectronaut_pic 1.5
#'@data upload :2024/08/16
#'

Spectronaut_pic <- function(dpi = 600,
                            width = 7,
                            height = 5.2){
  
  path<- file.choose()
  #返回上一级路径
  file.path <- dirname(path)
  setwd(file.path)
  # load packages
  message("loading packages...")
  pacman::p_load(UniProt.ws,tidyverse,readxl,data.table,patchwork,Peptides,ggrepel)
  # library(fuzzyjoin)
  rm(list=ls())
  # Modifications -----------------------------------------------------------
  modify <- data.table(modification = c("Lac","Pho","Pen","GG","Cr","Ac","Sc","3hb","SN","HNAc","Ph","Bu"),
                       mw = c(72.021129,-18.015287,112.01604399,114.0429274472,68.02621475,42.0105646,100.0160439947,86.03677943,107.004099,203.079372,79.96633,70.041864))
  # # 读取修饰文件夹
  # if (file.exists(modify.xlsx)) {  
  #   modify_add <- readxl::read_excel("modify.xlsx")
  #   modify <- rbind(modify,modify_add)
  # } else if (grepl("\\.csv$",filepath)) {  
  #   modify_add <- data.table::fread("modify.csv")
  #   modify <- rbind(modify,modify_add)
  # }
  
  # Input file --------------------------------------------------------------
  file <- dir() %>% as.data.table() %>% setnames(".","filename") %>% 
    .[(filename %like% "xls"|filename %like% "tsv") & !filename %like% "proteinlist"] %>% 
    .[,change:=apply(., 2, function(x){gsub("\\[Carbamidomethyl \\(C\\)\\]","",x)}) %>% as.data.table() %>% 
        # apply(., 2, function(x){gsub("\\[Oxidation \\(M\\)\\]","",x)}) %>% 
        apply(., 2, function(x){gsub(" \\(.+?\\)","",x)}) %>% 
        gsub("Y\\[Phospho\\]","Y\\[Ph\\]",.) %>% 
        gsub("M\\[Oxidation]","Z",.)] %>% 
    .[,`:=`(`modify sequence`=str_split(.$change,"_",simplify = T) %>% .[,2],
            modification = change %>% str_extract("\\[.+?\\]") %>% str_sub(2,-2))] %>%
    # .[,modification:=str_sub(modification,1,str_locate(modification," ") %>% .[,1]-1)] %>% 
    .[,sequence:=gsub("\\[.+?\\]","",.$`modify sequence`)] %>% 
    .[sequence!=""]
  file$modification <- gsub("Phospho","Pho",file$modification)
  file$modification <- gsub("Penta","Pen",file$modification)
  file$modification <- gsub("GlyGly","GG",file$modification)
  file$modification <- gsub("Crotonyl","Cr",file$modification)
  file$modification <- gsub("Acetyl","Ac",file$modification)
  file$modification <- gsub("Succinyl","Sc",file$modification)
  file$modification <- gsub("C8H12OS2","CS",file$modification)
  file$modification <- gsub("HexNAc","HNAc",file$modification)
  file$modification <- gsub("Butyrylation","Bu",file$modification)
  
  # fasta Analysis -------------------------------------------------------
  protein <- read_xlsx("proteinlist.xlsx") %>% .[,1] %>% as.data.table() %>% setnames(colnames(.),"Protein accession")
  # protein <- read_xls("proteinlist.xls") %>% .[,1] %>% as.data.table() %>% setnames(colnames(.),"Protein accession")
  
  sequence <- mapUniProt(
    from = "UniProtKB_AC-ID",
    to = "UniProtKB",
    columns = c("accession", "id", "gene_names","sequence"),
    query = list(ids = protein$`Protein accession`[!duplicated(protein$`Protein accession`)])
  ) %>% 
    dplyr::select(c("From","Entry.Name", "Gene.Names","Sequence")) %>% 
    setnames(colnames(.),c("uniprot","Entry.Name","gene","Sequence"))
  
  # Database ----------------------------------------------------------------
  AA <- aaList() %>% as.data.table() %>% setnames(colnames(.),"amino") %>% 
    .[,mw:=mw(amino,monoisotopic = T,aaShift = c(C= 57.021463))-18.01056] %>% 
    rbind(.,data.table(amino = "Z",
                       mw = (131.04049+15.994914)))
  
  dir.create("readtable")
  dir.create("Spectrogram")
  
  # Theoretical Mw ----------------------------------------------------------
  for (n in 1:nrow(file)) {
    theoretical <- file$`modify sequence`[n] %>% gsub("\\[.+?\\]","*",.) %>% str_split(pattern = "(?=[[:upper:]])",simplify = T) %>% .[-1] %>% as.data.table() %>% 
      setnames(".","amino") %>% 
      .[amino %like% "[*]", modification := file$modification[n]] %>% 
      .[,`:=`(b = paste("b",seq(1,nrow(.)),sep = ""),
              y = paste("y",seq(1,nrow(.)) %>% rev(),sep = ""))] %>% 
      .[(amino=="K"|amino=="R")&is.na(.$modification),charge:=1] %>%
      # .[(amino=="K"|amino=="R"|amino=="H")&is.na(.$modification),charge:=1] %>%
      .[1,charge:=1] %>% .[nrow(.),charge:=0] # b离子首尾电荷矫正
    theoretical$amino <- str_remove(theoretical$amino,"[*]")
    theoretical <- left_join(theoretical,AA)
    # theoretical$mw[1] <- 42.0105646+theoretical$mw[1]
    theoretical[modification==file$modification[n],mw:=mw+modify[modification==file$modification[n],mw]]
    theoretical$modification[is.na(theoretical$modification)] <- ""
    theoretical$charge[is.na(theoretical$charge)] <- 0
    # theoretical$modification[2] <- "Ox"
    # theoretical$mw[2] <- theoretical$mw[12]+15.9949146
    
    theoretical[,`:=`(charge=cumsum(theoretical$charge),
                      bmw=cumsum(theoretical$mw) + 18.015287 - 17.00735,
                      ymw=rev(theoretical$mw) %>% cumsum() %>% rev() + 18.015287 + 1.00794)]
    theoretical[b=="b1",`:=`(bmw=0,`m/z`=0)][y=="y1",`:=`(ymw=0,`m/z`=0)][b==paste("b",nrow(theoretical),sep = ""),`:=`(bmw=0,`m/z`=0)][y==paste("y",nrow(theoretical),sep = ""),`:=`(ymw=0,`m/z`=0)]
    
    yion <- theoretical %>% dplyr::select(c("amino","modification","y","charge","ymw")) %>% setnames(c("y","ymw"),c("ion","mw"))
    yion$charge <- max(yion$charge)-yion$charge+1
    
    theoretical <- rbind(theoretical %>% dplyr::select(c("amino","modification","b","charge","bmw")) %>% setnames(c("b","bmw"),c("ion","mw")),yion) %>% 
      .[,`:=`(`m/z`=(mw+(charge-1)*1.00794)/charge,
              type=str_sub(ion,1,1),
              label=paste(ion,str_dup("+",times = charge),sep = ""))]
    fragment <- theoretical %>% dplyr::select(c("label","m/z")) %>% 
      arrange(by = `m/z`) 
    if (grepl(".tsv",file[1,1])) {
      data <- fread(file[n,1] %>% as.character()) %>% 
        dplyr::select(ends_with("[Spectra]")) %>% 
        setnames(colnames(.),c("m/z","Intensity")) %>% 
        # .[Intensity > 10 & Intensity < 5000]
        # .[Intensity > 10 & Intensity < 5000][`m/z`<800]
        .[Intensity > 10]
    }else{
      file.copy(file$filename[n],"readtable/",overwrite = T)
      file.rename(paste("readtable/",file$filename[n],sep = ""),paste("readtable/",gsub(".xls",".csv",file$filename[n]),sep = ""))
      data <- fread(paste("readtable/",gsub(".xls",".csv",file$filename[n]),sep = "")) %>% 
        dplyr::select(ends_with("[Spectra]")) %>% 
        setnames(colnames(.),c("m/z","Intensity")) %>% 
        # .[Intensity > 10 & Intensity < 5000]
        # .[Intensity > 10 & Intensity < 5000][`m/z`<800]
        .[Intensity > 10]
    }
    
    
    # Match fragment ----------------------------------------------------------
    for (n1 in 1:nrow(fragment)) {
      data[abs(data$`m/z`-fragment[[n1,2]])<0.02,frag:=fragment$label[n1]]
    }
    data <- data[,group:= str_sub(data$frag,1,1) %>% as.factor()][!is.na(data$group),xlabel:= `m/z`][!is.na(data$group),ylabel:= frag %>% str_remove("[+]") %>% str_remove("[+]")] %>% 
      arrange(desc(Intensity))
    data <- rbind(data[!is.na(frag)] %>% unique(by="frag"),data[is.na(frag)])
    data$xlabel <- round(data$xlabel,2)
    
    # annotation --------------------------------------------------------------------
    annotation <- theoretical %>% dplyr::select(c("amino","modification")) %>% .[1:(nrow(theoretical)/2)] %>% 
      .[,`:=`(bion=paste("b",seq(1:nrow(.)),sep=""),
              yion=paste("y",seq(1:nrow(.)) %>% rev(),sep=""),
              x = seq(1,nrow(.)),y = 0)]
    annotation <- annotation[bion %in% data$ylabel,xlabel:=x][yion %in% data$ylabel,ylabel:=x]
    annotation <- annotation[amino=="Z",modification:="Ox"]
    annotation$amino <- gsub("Z","M",annotation$amino)
    
    # annotation$modification[1] <- "Ac"
    
    # Plot --------------------------------------------------------------------
    p1 <- ggplot(data,aes(x = `m/z`,y = Intensity))+
      geom_segment(aes(x = `m/z`, xend = `m/z`,y = 0,yend = Intensity),color = "aliceblue")+
      geom_segment(aes(x = `xlabel`, xend = `xlabel`,y = 0,yend = Intensity, color = group))+
      geom_text_repel(aes(label = frag, color = group),vjust = -0.4,hjust = 0.5,size = 3)+
      geom_text_repel(aes(label = xlabel, color = group,vjust = -1.5,hjust = 0.5),size = 3)+
      scale_color_manual(values = c("brown1","#0073C2FF")) +
      # scale_color_manual(values = c("brown1","#0073C2FF") %>% rev()) + 
      scale_y_continuous(expand = c(0,0),limits = c(0,data[group == "b"|group == "y"] %>% dplyr::select(.,"Intensity") %>% max(na.rm = T)*1.2))+
      ylab("Intensity")+
      theme_classic()+
      theme(legend.position = "none",aspect.ratio = 6/10)
    
    p2 <- ggplot(annotation,aes(x=x,y=y))+
      geom_text(aes(label = amino))+
      geom_text(aes(x=ylabel, y=y+0.25,label = yion),size = 3,color = "#0073C2FF")+
      geom_text(aes(x=xlabel, y=y-0.23,label = bion),size = 3,color = "brown1")+
      geom_text(aes(x=x, y=y+0.13,label = modification),size = 3, color = "grey2")+
      geom_segment(aes(x = xlabel-0.4, xend = xlabel+0.5, y =-0.3, yend = -0.3),color = "brown1")+
      geom_segment(aes(x = xlabel+0.5 , xend = xlabel+0.5, y = -0.3, yend = 0),color = "brown1")+
      geom_segment(aes(x = ylabel-0.5 , xend = ylabel-0.5, y = 0, yend = 0.3),color = "#0073C2FF")+
      geom_segment(aes(x = ylabel-0.5 , xend = ylabel+0.4, y = 0.3, yend = 0.3),color = "#0073C2FF")+
      coord_cartesian(ylim = c(-0.5,0.5))+
      theme_void()+
      theme(legend.position = "none",aspect.ratio = 1/5)
    
    p3 <- p2/p1
    
    # # path
    # if(!dir.exists(Results)){
    #   dir.create(Results)
    #   setwd(paste0(file.path,"/Results"))
    # }
    # output
    ggsave(paste("Spectrogram/",sequence$uniprot[file$sequence[n] %>% grep(sequence$Sequence)] %>% .[1],"_",file$`modify sequence`[n] %>% gsub("Z","M[Ox]",.),".pdf",sep = ""),
           p3,dpi = 600,width = 7,height = 5.2)
    ggsave(paste("Spectrogram/",sequence$uniprot[file$sequence[n] %>% grep(sequence$Sequence)] %>% .[1],"_",file$`modify sequence`[n] %>% gsub("Z","M[Ox]",.),".png",sep = ""),
           p3,dpi = 600,width = 7,height = 5.2)
    # fwrite(data,paste("Spectrogram/",sequensce$uniprot[file$sequence[n] %>% grep(sequence$Sequence)] %>% .[1],"_",file$`modify sequence`[n] %>% gsub("Z","M[Ox]",.),".csv",sep = ""))
  }
  unlink("readtable",recursive = T)
  message("DONE!!")
}

#'@title DIA_NN 
#'@data upload :2024/08/16
#'
DIANN_pic <- function(dpi = 600,
                      width = 7,
                      height = 5.2){
  
  path<- file.choose()
  #返回上一级路径
  file.path <- dirname(path)
  setwd(file.path)
  
  # load packages
  pacman::p_load(tidyverse,readxl,data.table,patchwork,ggrepel)
  rm(list=ls())
  sheets <- excel_sheets("pic.xlsx")
  lib <- fread("report-lib.tsv") %>% 
    unique(by = c("PeptideSequence","ProteinGroup")) %>% 
    dplyr::select(c("ProteinGroup","PeptideSequence"))
  filename <- read_xlsx("pic.xlsx",sheet = sheets[1],col_names = F) %>% .[1,1] %>% str_split("_",simplify = T) %>% .[1,1] %>% as.character()
  dir.create(filename)
  
  # Database ----------------------------------------------------------------
  AA <- data.table(amino = c("A","C","D","E","F","G","H","I","K","L","M","N","P","Q","R","S","T","V","W","Y","J","X","U","B","Z"),
                   mw = c(71.03711,103.00919+57.021463,115.02694,129.04259,147.06841,57.02146,137.05891,113.08406,128.09496,113.08406,
                          131.04049,114.04293,97.05276,128.05858,156.10111,87.03203,101.04768,99.06841,186.07931,163.06333,
                          173.03242,113.04768,113.08406,111.03203,103.00919))
  for (t in 1:length(sheets)) {
    # Data Input --------------------------------------------------------------
    data <- read_xlsx("pic.xlsx",skip = 2,sheet = sheets[t]) %>% as.data.table() %>% 
      .[,Intensity := Intensity/10000] %>% 
      .[Intensity>0.00001]        # 调整阈值
    # .[Intensity>0.0001&Intensity<100]
    name <- read_xlsx("pic.xlsx",sheet = sheets[t]) %>% .[1,1] %>% str_split(",",simplify = T) %>% .[1,1]
    
    # Theoretical Mw ----------------------------------------------------------
    theoretical <- str_split(name,pattern = "",simplify = T) %>% t() %>% as.data.table() %>% 
      setnames("V1","amino") %>% 
      .[,`:=`(b = paste("b",seq(1,nrow(.)),sep = ""),
              y = paste("y",seq(1,nrow(.)) %>% rev(),sep = ""))] %>% 
      .[(amino=="K"|amino=="R")&is.na(.$modification),charge:=1] %>% 
      .[1,charge:=1] %>% .[nrow(.),charge:=0] %>% # b离子首尾电荷矫正
      left_join(.,AA)
    
    # theoretical$mw[site] <- theoretical$mw[site]-18.015287
    
    theoretical$charge[is.na(theoretical$charge)] <- 0
    theoretical[,`:=`(charge=cumsum(theoretical$charge),
                      bmw=cumsum(theoretical$mw) + 18.015287 - 17.00735,
                      ymw=rev(theoretical$mw) %>% cumsum() %>% rev() + 18.015287 + 1.00794)]
    theoretical[b=="b1",`:=`(bmw=0,`m/z`=0)][y=="y1",`:=`(ymw=0,`m/z`=0)][
      b==paste("b",nrow(theoretical),sep = ""),`:=`(bmw=0,`m/z`=0)][
        y==paste("y",nrow(theoretical),sep = ""),`:=`(ymw=0,`m/z`=0)]
    
    yion <- theoretical %>% dplyr::select(c("amino","y","charge","ymw")) %>% setnames(c("y","ymw"),c("ion","mw"))
    yion$charge <- max(yion$charge)-yion$charge+1
    
    theoretical <- rbind(theoretical %>% dplyr::select(c("amino","b","charge","bmw")) %>% setnames(c("b","bmw"),c("ion","mw")),
                         theoretical %>% dplyr::select(c("amino","y","charge","ymw")) %>% setnames(c("y","ymw"),c("ion","mw"))) %>% 
      .[,`:=`(`m/z`=(mw+(charge-1)*1.00794)/charge,
              type=str_sub(ion,1,1),
              frag=paste(ion,str_dup("+",times = charge),sep = ""))]
    fragment <- theoretical %>% dplyr::select(c("frag","m/z")) %>% arrange(by = `m/z`)
    
    # Match fragment ----------------------------------------------------------
    for (n in 1:nrow(fragment)) {
      data <- data[abs(data$`m/z`-fragment[[n,2]]) < 0.02,frag:= fragment[n,1]]   #调整精度
    }
    data <- data[,group:= str_sub(data$frag,1,1) %>% as.factor()][!is.na(data$group),xlabel:= `m/z`][!is.na(data$group),ylabel:= frag %>% str_remove("[+]") %>% str_remove("[+]")] %>% 
      arrange(desc(Intensity))
    data <- rbind(data[!is.na(frag)] %>% unique(by="frag"),data[is.na(frag)])
    data$xlabel <- round(data$xlabel,2)
    sup_data <- data.table(`m/z`=c(200,200),Intensity=c(0,0),frag=NA,group=c("b","y"),xlabel=NA,ylabel=NA)
    data <- rbind(data,sup_data)
    
    # annotation --------------------------------------------------------------------
    annotation <- theoretical %>% dplyr::select(c("amino")) %>% .[1:(nrow(theoretical)/2)] %>% 
      .[,`:=`(bion=paste("b",seq(1:nrow(.)),sep=""),
              yion=paste("y",seq(1:nrow(.)) %>% rev(),sep=""),
              x = seq(1,nrow(.)),y = 0)]
    annotation <- annotation[bion %in% data$ylabel,xlabel:=x][yion %in% data$ylabel,ylabel:=x]
    # annotation <- annotation[,modification:= ""]
    # annotation$modification[site] <- "pho"
    
    # Plot --------------------------------------------------------------------
    p1 <- ggplot(data,aes(x = `m/z`,y = Intensity))+
      geom_segment(aes(x = `m/z`, xend = `m/z`,y = 0,yend = Intensity),color = "aliceblue")+
      geom_segment(aes(x = xlabel, xend = xlabel,y = 0,yend = Intensity, color = group))+
      geom_text(aes(label = frag, color = group),vjust = -0.4,hjust = 0.5,size = 3)+
      geom_text(aes(label = xlabel, color = group,vjust = -1.5,hjust = 0.5),size = 3)+
      scale_color_manual(values = c("#0073C2FF","brown1") %>% rev()) + 
      # scale_y_continuous(expand = c(0,0),limits = c(0,max(data$Intensity)*1.2))+
      scale_y_continuous(expand = c(0,0),limits = c(0,data[group=="b"|group=="y"] %>% dplyr::select(Intensity) %>% max(na.rm = T)*1.5))+
      ylab("Intensity [counts] (10^5)")+
      theme_classic()+
      theme(legend.position = "none",aspect.ratio = 6/10)
    # 
    p2 <- ggplot(annotation,aes(x=x,y=y))+
      geom_text(aes(label = amino))+
      geom_text_repel(aes(x=ylabel, y=y+0.25,label = yion),size = 3,color = "#0073C2FF")+
      geom_text_repel(aes(x=xlabel, y=y-0.23,label = bion),size = 3,color = "brown1")+
      geom_segment(aes(x = xlabel-0.4, xend = xlabel+0.5, y =-0.3, yend = -0.3),color = "brown1")+
      geom_segment(aes(x = xlabel+0.5 , xend = xlabel+0.5, y = -0.3, yend = 0),color = "brown1")+
      geom_segment(aes(x = ylabel-0.5 , xend = ylabel-0.5, y = 0, yend = 0.3),color = "#0073C2FF")+
      geom_segment(aes(x = ylabel-0.5 , xend = ylabel+0.4, y = 0.3, yend = 0.3),color = "#0073C2FF")+
      coord_cartesian(ylim = c(-0.5,0.5))+
      theme_void()+
      theme(legend.position = "none",aspect.ratio = 1/5)
    
    p3 <- p2/p1
    
    # # path
    # if(!dir.exists(Results)){
    #   dir.create(Results)
    #   setwd(paste0(file.path,"/Results"))
    # }
    # output
    ggsave(paste(filename,"/",lib[PeptideSequence==name] %>% .[1,1] %>% as.character(),"_",name,".pdf",sep = ""),p3,dpi = 600,width = 7,height = 5.2)
    ggsave(paste(filename,"/",lib[PeptideSequence==name] %>% .[1,1] %>% as.character(),"_",name,".png",sep = ""),p3,dpi = 600,width = 7,height = 5.2)
    
  }
  message("DONE")
}
#'@title DIA_Fig 
#'@data upload :2024/08/16
#'

DIA_Fig <- function(dpi = 600,
                    width = 7,
                    height = 5.2){
  path<- file.choose()
  #返回上一级路径
  file.path <- dirname(path)
  setwd(file.path)
  
  if (file.exists("report-lib.tsv")) {
    message("frome DIA-NN...strating")
    DIANN_pic(dpi = dpi,width = width,height = height)
  }else{
    message("frome Spectronaut...strating")
    Spectronaut_pic(dpi = dpi,width = width,height = height)
  }
  message("DONE!")
}

#'@title spectronaut qc

spectronaut_PTM_depth <- function(filepath = NULL,
                                  rm_modify = NULL,
                                  SiteProbability = 0.75){
  
  if(is.null(filepath)){
    path <- file.choose()
  }else{
    path <- filepath
  }
  #检查文件选择是否正确
  if(!path %like% "-PTMSiteReport.tsv"){
    stop("文件选择错误!")
  }
  dat <- data.table::fread(path)
  
  #全部
  #筛选修饰 和位点准确性
  if(is.null(rm_modify)){
    rm_mod <- "NA"
  }else{
    rm_mod <- paste(rm_modify,collapse = "|")
  }
  
  df <- dat %>% filter(!grepl(rm_mod, PTM.SiteAA, ignore.case = TRUE)) %>% 
    filter(PTM.SiteProbability >= SiteProbability)
  #修饰类型
  modify <- unique(df$PTM.ModificationTitle)

  Results <- data.frame()
  
  for (i in modify) {
    
    df1 <- df %>% filter(PTM.ModificationTitle == i)
    
    #unique 总深度
    df2 <- data.frame(df1$PG.ProteinNames,df1$PTM.SiteLocation) %>% unique()
    
    result1 <- data.frame(info = "all", 
                          value = nrow(df2),
                          modify = i)
    # 样本类型
    samples <- unique(dat$R.Condition)
    
    if("Not Defined" %in% samples){
      
      Results <- rbind(Results,result1)
      
    }else{
      
      result2 <- lapply(samples, function(x){
        
        #筛选修饰 和大于0.75
        sample_df <- dat %>% filter(R.Condition == x) %>% 
          filter(!grepl(rm_mod, PTM.SiteAA, ignore.case = TRUE)) %>% 
          filter(PTM.SiteProbability >= SiteProbability)
        
        sample_df <- sample_df %>% filter(PTM.ModificationTitle == i)
        
        #unique 总深度
        dat2 <- data.frame(sample_df$PG.ProteinNames,sample_df$PTM.SiteLocation) %>% unique()
        
        data.frame(info = x, 
                   value = nrow(dat2),
                   modify = i)
        
      }) %>% do.call(rbind,.)
      Results <- rbind(Results,result1,result2)
    }
  }
  
  return(Results)
  
}

#########################################
#需要先设置默认路径

mult_spectronaut_PTM_depth <- function(write = FALSE){
  
  tsv.paths <- dir(pattern = "\\-PTMSiteReport.tsv",recursive = TRUE,full.names = TRUE)
  
  Results <- data.frame()
  
  for (i in tsv.paths) {
    
    df <- ProteomicsR::spectronaut_PTM_depth(filepath = i)
    
    message("\033[32m")
    print(i)
    print(df)
    message("-----------------------------------")
    
    dfRes <- cbind(df,tsv.path=i)
    
    Results <- rbind(Results,dfRes)
    
    if(write){
      data.table::fwrite(Results,"depth.csv")
    }
  }
}



