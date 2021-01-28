#Micheal Taylor microarray set Investigation

library(limma)
library(oligo)
library(siggenes)
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(tibble)




######################################################################################################
#Loading CEL files into R from local machine for Medulloblastoma Patients#############################-
######################################################################################################
setwd("C:/Users/12298/Desktop/Data_Analytics/Taylor_2017/CEL_Files/GSE85217_RAW")
#CELfiles=dir(pattern="CEL.gz$")
#rawdata_MB=read.celfiles(filenames= CELfiles[1:763])

#Normalizing and Summarizing expression data in Annotated Expression object/Extraction from object using exprs function
#norm_MB_dat=rma(rawdata_MB)
#Expr_norm=exprs(norm_MB_dat)
#saveRDS(Expr_norm, "Normalized_MB.rds") #saving normalized file for use later on
Expr_norm <- readRDS("Normalized_MB.rds")
colnames(Expr_norm)=substr(colnames(Expr_norm), 1,10)

#Boxplot following normalization to inspect visually
boxplot(Expr_norm[1:nrow(Expr_norm),1:50])

#Loading Affy Probe IDs - Bioconductor package was missing probes - used txt supplied with CEL files
setwd("C:/Users/12298/Desktop/Data_Analytics/Taylor_2017")
Probes_Final <- read.table("GPL6244-17930.txt", fill=TRUE)
Probes_Final <- Probes_Final[,c(1,2,12)] %>% janitor::row_to_names(row_number = 1) %>% `colnames<-`(c("ID", "RefSeq", "Gene_Name"))


#Bringing in patient info and cleaning up
setwd("C:/Users/12298/Desktop/Data_Analytics/Taylor_2017")
Patient_Info=read.csv("GSE85219_series_matrix.csv")
Patient_Info=Patient_Info[c(32, 41, 42),] %>% t() %>% `colnames<-`(c("Patient_ID", "Subgroup", "Subtype"))
rownames(Patient_Info) <- NULL
Patient_Info=Patient_Info[-1,]


####################################################################################################
#Reading in CEL files for cerebellar controls and normalizing#######################################-
####################################################################################################
setwd("c:/Users/12298/Desktop/Data_Analytics/Normal_Cerebellum/GSE7307_RAW/")
CELfiles=dir(pattern="CEL.gz")
rawdata_CB=read.celfiles(filenames= CELfiles[1:8])

#Performing RMA normalization and then extracting log2 expression values 
norm_CB_dat = rma(rawdata_CB)
Expr_norm=exprs(norm_CB_dat)

saveRDS(Expr_norm, "Normalized_CB.rds") #saving normalized file for use later on
Expr_norm <- readRDS("Normalized_CB.rds")
colnames(Expr_norm)=substr(colnames(Expr_norm), 1,10)

#Inspecting normalization
boxplot(Expr_norm[1:nrow(Expr_norm),1:8])


#Attempting to pull probes from bioconductor for normal cerebellum array - must call names on input matrix - format requires probes as colnames
library(hgu133a.db)

Expr_norm <- t(Expr_norm)
commonGenes=colnames(Expr_norm)
geneNames=as.character(hgu133aACCNUM[commonGenes])
head(geneNames)

setwd("C:/Users/12298/Desktop/Data_Analytics/Normal_Cerebellum/")
Probes_Final <- read.table("GPL570-55999.txt", fill=TRUE)
Probes_Final <- Probes_Final[,c(1,2,12)] %>% janitor::row_to_names(row_number = 1) %>% `colnames<-`(c("ID", "RefSeq", "Gene_Name"))



#Filter based on common gene names between control cerebellum and MB patients




#Extract Housekeeping genes and compare between CB and MB averages


##################################################################################################
#Creating YBX plots and comparing with normal cerebellum##########################################-
##################################################################################################
#Expr_norm <- data.table::setDT(as.data.frame(Expr_norm), keep.rownames=TRUE)[]
#colnames(Expr_norm)[1] <- "ID"
#Merge_DF <- left_join(as.data.frame(Probes_Final), as.data.frame(Expr_norm), by="ID", incomparables=NULL)
#saveRDS(Merge_DF, "Finalized_Expression_Array.rds")
Merge_DF <- readRDS("Finalized_Expression_Array.rds")

#Extracting relevant genes and creating data frame
YBX1_Correlation <- c("YBX1", "YBX2", "YBX3")

YBX_Panel <- Merge_DF[Merge_DF$Gene_Name %in% YBX1_Correlation, -(1:2)] %>% pivot_longer(values_to="2Log_Expression", names_to = "Patient_ID", -Gene_Name)  
YBX_Graph <- c("YBX1", "YBX3", "YBX2")
YBX_Panel <- Merge_DF[Merge_DF$Gene_Name %in% YBX_Graph, -(1:2)] %>% pivot_longer(values_to="2Log_Expression", names_to = "Patient_ID", -Gene_Name)  

#Merging patient data with genes extracted
YBX_Panel <- merge(YBX_Panel, Patient_Info, by="Patient_ID")
YBX_Panel$Subgroup <- gsub("^.{0,10}", "", YBX_Panel$Subgroup)
YBX_Panel$Subtype <- gsub("^.{0,9}", "", YBX_Panel$Subtype)


#Geom_Point Plots for YBX Probes 
ggplot(YBX_Panel, aes(x=Subgroup, y=`2Log_Expression`, color=Gene_Name)) + geom_boxplot() + labs(title="YBX Expression Across Medulloblastoma Subgroups") + ylab("2Log Expression") + ggsave("YBX_in_MB_Subgroups.png")

YBX_Panel %>% filter(Subgroup == "SHH") %>% ggplot(aes(x=Subtype, y=`2Log_Expression`, color=Gene_Name)) + geom_point() + geom_jitter(width=0.35) + stat_summary(fun=mean, geom="line", color="black", aes(group=Gene_Name)) + labs(title="YBX Expression in SHH Subtypes") + ggsave("YBX_in_SHH_Subtypes.png")



##################################################################################################
#Creating list of correlation coefficients by looping#############################################-
##################################################################################################
#Averaging between duplicate probes, pivoting longer, merging with patient ID, filtering for SHH patients only, filtering for 'intact' values, and pivoting wider before looping
Merge_DF <- readRDS("Finalized_Expression_Array.rds")
YBX_Panel <- Merge_DF[,-c(1:2)]

YBX_Panel <- YBX_Panel[!(YBX_Panel$Gene_Name == '' | YBX_Panel$Gene_Name == '//' | YBX_Panel$Gene_Name == '//' | YBX_Panel$Gene_Name == 'NONCODE' | YBX_Panel$Gene_Name == '---' | YBX_Panel$Gene_Name == '///' | YBX_Panel$Gene_Name == 'ENSEMBL'),]



YBX_Panel <- aggregate(YBX_Panel[-1], YBX_Panel[1], mean) #Averaging Duplicate probes
YBX_Panel <- YBX_Panel %>% pivot_longer(values_to="Log_Expression", names_to = "Patient_ID", -Gene_Name)
YBX_Panel <- YBX_Panel[!(is.na(YBX_Panel$Log_Expression)),] #Removing NA values from Expression Column

YBX_Panel <- merge(YBX_Panel, Patient_Info, by="Patient_ID")
YBX_Panel$Subgroup <- gsub("^.{0,10}", "", YBX_Panel$Subgroup)
YBX_Panel$Subtype <- gsub("^.{0,9}", "", YBX_Panel$Subtype)
YBX_Panel <- YBX_Panel %>% filter(Subgroup == "SHH")
saveRDS(YBX_Panel, "YBX_Long.rds")
YBX_Long <- readRDS("YBX_Long.rds")

YBX_Cor_Test <-  pivot_wider(YBX_Panel, names_from = "Patient_ID", values_from = "Log_Expression", -c("Subgroup", "Subtype"))
saveRDS(YBX_Cor_Test, "YBX_DF_For_Cor_Test")
YBX_Cor_Test <- readRDS("YBX_DF_For_Cor_Test")

YBX_Cor_Test_2 <- YBX_Cor_Test %>% remove_rownames %>% tibble::column_to_rownames(var="Gene_Name")

#Creating empty vector to store output
Cor_Vector_YBX1 <- c()

  #Looping
for (i in 1:nrow(YBX_Cor_Test_2))
  {
  print(i)
  Cor_Vector_YBX1[i] = cor(as.numeric(YBX_Cor_Test_2[17970,]), as.numeric   (YBX_Cor_Test_2[i,]))
  }

YBX_Cor_Results <- cbind(YBX_Cor_Test, Cor_Vector_YBX1)
YBX_Cor_Results <- as.data.frame(YBX_Cor_Results[,-c(2:224)])
YBX_Cor_Results <- rownames_to_column(YBX_Cor_Results, var = "Gene_Name")


saveRDS(YBX_Cor_Results, "YBX1_Cor_Results")
YBX_Cor_Results <- readRDS("YBX1_Cor_Results")
setwd("C:/Users/12298/Desktop/Data_Analytics/Taylor_2017")

cor(as.numeric(YBX_Cor_Test_2[17970,]), as.numeric(YBX_Cor_Test_2[16916,]))
plot(as.numeric(YBX_Cor_Test_2[17970,]), as.numeric(YBX_Cor_Test_2[16916,]), ylab="PARP1", xlab="YBX1", )

YBX1 <- YBX_Cor_Test_2[c(17970,16916),] %>% t() %>% as.data.frame()

ggplot(YBX1, aes(x=YBX1,y=PARP1)) + geom_point(pch=16, size=2, alpha=1, color = "#EF70DE") + ggtitle("YBX1-PARP1 Correlation in MB SHH Patients") + geom_smooth(method=lm, color="skyblue", fill="black", alpha=0.5) + geom_text(x=10.6, y=10.9, label="r2 = 0.6366", size=5) + theme_minimal(base_size = 16) + theme( legend.position = "none")

#End-----





