#Load the R Packages 


library("matrixStats") # version 0.58.0

library("data.table") # version 1.14.3

library(dplyr)

#Initialize variable

#Tissue name (Add all organs you want to compare)
TissTitle=c("Cecum","Duodenum","Distal_Colon", "Esophagus", "Heart_A", "Heart_B", "Heart_C", "Heart_D", "Ileum", "Jejunum", "Proximal_Colon", "Upper_Stomach")

#Comparison list (FC between two consecutive endpoints)
comg=c("newratio_3v0","newratio_12v3")

main_dir <- "~/MElaS11.21" #search for main directory

nested_dir <- list.dirs(main_dir, recursive = FALSE) #search nested directories

for (path in nested_dir) {
  print(paste("Processing directory:", main_dir))
  setwd(path)
  
  #Input the expression dynamics and statistical significance information
  for(tissi in TissTitle){
    for(tcomgi in comg){
      tFCFDR=fread(paste0(tissi,"/",tcomgi,".txt"),sep = "\t",data.table = F)
      rownames(tFCFDR)=tFCFDR$ID
      assign(paste0(tissi,".",tcomgi,".FCFDR"),tFCFDR)
    }
  }
  #Calculate the Metabolic Elasticity Score (MElaS)
  for(tissi in TissTitle){
    #Variable name for Comparison
    tFN.name=paste0(tissi,".newratio_3v0.FCFDR")
    tRF.name=paste0(tissi,".newratio_12v3.FCFDR")
    
    #Create the data.frame to store the MElaS
    tScoreMat=data.frame(ID=get(tFN.name)[,"ID"])
    
    #Metabolite List
    tgid=tScoreMat$ID
    
    #Expression dynamics 
    tFN.FC = get(tFN.name)[tgid,"FC"]
    tRF.FC = get(tRF.name)[tgid,"FC"]
    
    #Sign of expression dynamics
    tMultilog2FCSign = -1*sign(tFN.FC * tRF.FC)
    tMultilog2FCSign[tMultilog2FCSign <= 0] = 0
    
    #Absolute expression dynamics matrix
    tAbslog2FC.Mat=data.frame(Abslog2FC.FN=abs(get(tFN.name)[tgid,"FC"]),
                              Abslog2FC.RF=abs(get(tRF.name)[tgid,"FC"]),
                              stringsAsFactors = F)
    
    #Calculate restoration extent of gene expression
    tlog2FC.Ratio=rowMins(as.matrix(tAbslog2FC.Mat))/rowMaxs(as.matrix(tAbslog2FC.Mat))
    tlog2FC.Ratio[is.na(tlog2FC.Ratio)]=0
    
    #Calculate the weight for FDR (statistical significance)
    #tFDR.FN=log10(get(tFN.name)[tgid,"adj.P.Val"])
    tFDR.RF=log10(get(tRF.name)[tgid,"adj.P.Val"])
    #tFDR.FN[tFDR.FN >= log10(0.05)]=tFDR.FN[tFDR.FN >= log10(0.05)]/(log10(0.05))
    #tFDR.FN[tFDR.FN < log10(0.05)]=1
    tFDR.RF[tFDR.RF >= log10(0.05)]=tFDR.RF[tFDR.RF >= log10(0.05)]/(log10(0.05))
    tFDR.RF[tFDR.RF < log10(0.05)]=1
    
    #Integrate the expression dynamics, statistical significance, and restoration extent to 
    #calculate the MElaS
    tAbslog2FC.FDRSum=(tAbslog2FC.Mat$Abslog2FC.FN+tAbslog2FC.Mat$Abslog2FC.RF*tFDR.RF)
    tGElaS= tMultilog2FCSign*tlog2FC.Ratio*tAbslog2FC.FDRSum
    tScoreMat[,paste0("MElaS")]=data.frame(tGElaS,stringsAsFactors = F)
    
    #Output the MElaS
    fwrite(tScoreMat,file = paste0(tissi,"/", tissi, "_MElaS.txt"),sep = "\t")
    
    
  }
  #Add LogFC of 3v0 and 12v3 to MElaS file
  for(tissi in TissTitle){
    melas_df <- read.table(paste0(tissi, "/", tissi, "_MElaS.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    logfc_3v0_df <- read.table(paste0(tissi, "/", "newratio_3v0.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    logfc_12v3_df <- read.table(paste0(tissi, "/", "newratio_12v3.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    
    
    combined_df <- melas_df %>%
      left_join(logfc_3v0_df %>% select(ID, logFC_3v0=FC), by = "ID") %>%
      left_join(logfc_12v3_df %>% select(ID, logFC_12v3=FC), by = "ID")
    
    
    write.table(combined_df, paste0(tissi, "/", tissi, "_MElaS.txt"), sep = "\t", row.names = FALSE, quote = FALSE)
    
    
    #Add LogFC of 0w, 3w and 12w to MElaS file
    melas_df <- read.table(paste0(tissi, "/", tissi, "_MElaS.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
    logfc_0w_df <- read.csv(paste0(tissi, "/", "newratio_0w.csv"), header = TRUE, stringsAsFactors = FALSE)
    logfc_3w_df <- read.csv(paste0(tissi, "/", "newratio_3w.csv"), header = TRUE, stringsAsFactors = FALSE)
    logfc_12w_df <- read.csv(paste0(tissi, "/", "newratio_12w.csv"), header = TRUE, stringsAsFactors = FALSE)
    logfc_24w_df <- read.csv(paste0(tissi, "/", "newratio_24w.csv"), header = TRUE, stringsAsFactors = FALSE)
    
    
    combined_df <- melas_df %>%
      left_join(logfc_0w_df %>% select(ID, FC_0w=FC), by = "ID") %>%
      left_join(logfc_3w_df %>% select(ID, FC_3w=FC), by = "ID") %>%
      left_join(logfc_12w_df %>% select(ID, FC_12w=FC), by = "ID") %>%
      left_join(logfc_24w_df %>% select(ID, FC_24w=FC), by = "ID")
    
    
    write.table(combined_df, paste0(tissi, "/", tissi, "_MElaS.txt"), sep = "\t", row.names = FALSE, quote = FALSE)

  }
  
}

#combines all MElaS dataframes into one file that includes all parasite species and infected parasites
datalist = list()
for (path in nested_dir) {
  print(paste("Processing directory:", path))
  for (organ in TissTitle) {
    print(organ)
    file_path <- file.path(path, organ, paste0(organ, "_MElaS.txt"))
    df <- read.table(file_path, header = TRUE, sep = "\t")
    
    df$organ <- organ
    #df$parasite_species<-substr(path, start = 32, stop = 33)
    df$parasite_species <- str_extract(path, "(?<=/MElaS11.21/MElaS).*?(?=0_100)")
    df$infected_parasites<-sub(".*_", "", path)
    #df$time<-time
    datalist[[length(datalist) + 1]] <- df
  }
}

combined_df<-do.call(rbind, datalist)
write.csv(combined_df, "All_DATA_2.3.csv", row.names = FALSE)
