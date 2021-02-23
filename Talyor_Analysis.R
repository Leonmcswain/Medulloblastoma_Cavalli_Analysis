#Cavalli microarray set Investigation accession GSE85217 paired with normal cerebellum from GSE44971

#Libraries for Microarray
library(limma)
library(oligo)
library(siggenes)

#Libraries for Data Munging and graphing 
library(dplyr)
library(ggplot2)
library(tidyr)
library(data.table)
library(tibble)


####################################################################################################
####################----------Medulloblastoma Patient Data --------------###########################-
####################################################################################################
setwd("C:/Users/12298/Desktop/Data_Analytics/Taylor_2017/CEL_Files/GSE85217_RAW")
#CELfiles=dir(pattern="CEL.gz$")
#rawdata_MB=read.celfiles(filenames= CELfiles[1:763])

#Normalizing and Summarizing expression data in Annotated Expression object/Extraction from object using exprs function
norm_MB_dat=rma(rawdata_MB)
Expr_norm=exprs(norm_MB_dat)
saveRDS(Expr_norm, "Normalized_MB.rds") #saving normalized file for use later on
Expr_norm_MB <- readRDS("C:/Users/12298/Desktop/Data_Analytics/Taylor_2017/CEL_Files/GSE85217_RAW/Normalized_MB.rds")
colnames(Expr_norm_MB)=substr(colnames(Expr_norm_MB), 1,10)

#Boxplot following normalization to inspect visually
boxplot(Expr_norm_MB[1:nrow(Expr_norm_MB),1:50])

#Bringing probes from txt file and merging with expression data 
Probes_Final <- read.table("C:/Users/12298/Desktop/Data_Analytics/Taylor_2017/GPL6244-17930.txt", fill=TRUE)
Probes_Final <- Probes_Final[,c(1,12)] %>% janitor::row_to_names(row_number = 1) %>% `colnames<-`(c("ID", "Gene_Name")) %>% as.data.frame()
Expr_norm_MB <- data.table::setDT(as.data.frame(Expr_norm_MB), keep.rownames=TRUE)[]
colnames(Expr_norm_MB)[1] <- "ID"
Merge_DF <- left_join(Probes_Final, Expr_norm_MB, by="ID", incomparables=NULL)
saveRDS(Merge_DF, "C:/Users/12298/Desktop/Data_Analytics/Taylor_2017/Finalized_Expression_Array.rds")


####################################################################################################
###################################------Cerebellar Controls----------##############################-
####################################################################################################
#GEO_Accession: GSE44971 - YBX1 values are not consistent here while genes in other rows are. Also, average is around 11 which is no different than tumor tissue. This is reminiscent of qPCR data I have gotten in the lab showing no differential transcript levels but differential protein levels. 
setwd("C:/Users/12298/Desktop/Data_Analytics/Normal_Cerebellum/GSE22569_RAW/")
CELfiles=dir(pattern="CEL.gz")
rawdata_CB=read.celfiles(filenames= CELfiles[1:10])

#Performing RMA normalization and then extracting log2 expression values 
norm_CB_dat = rma(rawdata_CB)
Expr_norm_CB <- exprs(norm_CB_dat)

#Inspecting normalization
boxplot(Expr_norm_CB[1:nrow(Expr_norm_CB),1:8])

#Merging probes with Expression data using accession number 
Probes_Final <- read.table("C:/Users/12298/Desktop/Data_Analytics/Normal_Cerebellum/GSE22569_RAW/GPL6244-17930.txt", fill=TRUE)
Probes_Final <- Probes_Final[,c(1,12)] %>% janitor::row_to_names(row_number = 1) %>% `colnames<-`(c("ID", "Gene_Name"))
Expr_norm_CB <- data.table::setDT(as.data.frame(Expr_norm_CB), keep.rownames=TRUE)[]
colnames(Expr_norm_CB)[1] <- "ID"
Merge_DF <- left_join(Probes_Final, Expr_norm_CB, by="ID", incomparables=NULL)
saveRDS(Merge_DF, "C:/Users/12298/Desktop/Data_Analytics/Normal_Cerebellum/GSE22569_RAW/Finalized_Expression_Array.rds")



####################################################################################################
#########-------Merging and Checking Housekeeping Genes----------###################################-
####################################################################################################
#Merge Expr_norm_CB and Expr_norm_MB by probe
Expr_norm_MB <- readRDS("C:/Users/12298/Desktop/Data_Analytics/Taylor_2017/Finalized_Expression_Array.rds") 
Expr_norm_CB <- readRDS("C:/Users/12298/Desktop/Data_Analytics/Normal_Cerebellum/GSE22569_RAW/Normalized_CB.rds") %>% as.data.frame()
Expr_norm_CB$ID <- rownames(Expr_norm_CB)
Merge_DF <- merge(Expr_norm_CB, Expr_norm_MB, by="ID")
saveRDS(Merge_DF, "C:/Users/12298/Desktop/Data_Analytics/Taylor_2017/Merge_Normal_Cerebellum.rds")

#Filter based on common gene names between control cerebellum and MB patients - E. Eisenberg and E.Y. Levanon, Trends in Genetics, 29 (2013) for list of housekeeping genes (~3000)
Housekeep_Genes <- rio::import(file="https://www.tau.ac.il/~elieis/HKG/HK_exons.xlsx") 
Housekeep_Genes <- unique(Housekeep_Genes$"Gene Name")
Housekeep_Gin <- filter(Merge_DF, Merge_DF$Gene_Name %in% Housekeep_Genes)
plot(rowMeans(Housekeep_Gin[,2:11]), rowMeans(Housekeep_Gin[,13:ncol(Housekeep_Gin)]))
cor(rowMeans(Housekeep_Gin[,2:11]), rowMeans(Housekeep_Gin[,13:ncol(Housekeep_Gin)]))

#Correlation is moderately high. I will deem this as acceptable. With >3000 genes I am sure there are bound to be deferentially expressed genes in there. May choose a smaller list of more acceptable housekeeping genes. 


####################################################################################################
#########-------Creating data frame with Normal Cerebellum and MB Patients--------##################-
####################################################################################################
Merge_DF <- readRDS("C:/Users/12298/Desktop/Data_Analytics/Taylor_2017/Merge_Normal_Cerebellum.rds")
#Creating data frame for Normal Cerebellum IDs to merge with Patient IDs before merging with expression values l
Normal_CB_DF <- data.frame("Patient_ID"=colnames(Expr_norm_CB[,1:10]), "Subgroup"=rep("Normal", 10), "Subtype"=rep("Normal",10))

#Bringing in patient info and cleaning up - bind with normal cerebellum mock data frame
Patient_Info=read.csv("C:/Users/12298/Desktop/Data_Analytics/Taylor_2017/GSE85219_series_matrix.csv")
Patient_Info=Patient_Info[c(32, 41, 42),] %>% t() %>% `colnames<-`(c("Patient_ID", "Subgroup", "Subtype"))
rownames(Patient_Info) <- NULL
Patient_Info <- Patient_Info[-1,] %>% as.data.frame()
Patient_Info$Subgroup <- gsub("^.{0,10}", "", Patient_Info$Subgroup)
Patient_Info$Subtype <- gsub("^.{0,9}", "", Patient_Info$Subtype)
Patient_Info <- rbind(Normal_CB_DF, Patient_Info)

#Removing rows with no gene identifier - microarray controls and non-coding genes
Merge_DF <- Merge_DF[!(Merge_DF$Gene_Name == '' | Merge_DF$Gene_Name == '//' | Merge_DF$Gene_Name == '//' | Merge_DF$Gene_Name == 'NONCODE' | Merge_DF$Gene_Name == '---' | Merge_DF$Gene_Name == '///' | Merge_DF$Gene_Name == 'ENSEMBL'),]

#Pivot Expression data frame
Long_Expression <- Merge_DF %>% pivot_longer(values_to="2Log_Expression", names_to = "Patient_ID", -c(Gene_Name, ID))

#Merging patient data with genes extracted
Gene_Extraction <- merge(Long_Expression, Patient_Info, by="Patient_ID")
saveRDS(Gene_Extraction, "C:/Users/12298/Desktop/Data_Analytics/Taylor_2017/CB_MB_Final.rds")

#################################################################################################
###################### Gene Extraction and Plotting #############################################-
#################################################################################################
Gene_Extraction <- readRDS("C:/Users/12298/Desktop/Data_Analytics/Taylor_2017/CB_MB_Final.rds")


#Geom_Point Plots for YBX Probes 

arrange(Gene_Extraction) %>% 
  mutate(Subgroup= factor(Subgroup, levels=c("Normal", "Group3", "Group4", "SHH", "WNT"))) %>%
ggplot(aes(x=Subgroup, y=`2Log_Expression`, fill=Gene_Name)) + geom_boxplot() + labs(title="YBX Expression Across Medulloblastoma Subgroups") + ylab("2Log Expression") + theme_classic(base_size=15) + ggsave("YBX_in_MB_Subgroups.png")

Gene_Extraction %>% filter(Subgroup == "SHH") %>% ggplot(aes(x=Subtype, y=`2Log_Expression`, color=Gene_Name)) + geom_point() + geom_jitter(width=0.35) + stat_summary(fun=mean, geom="line", color="black", aes(group=Gene_Name)) + labs(title="YBX Expression in SHH Subtypes") + ggsave("YBX_in_SHH_Subtypes.png")



####################################################################################################
#########-------Creating list of correlation coefficients by looping-----------#####################-
####################################################################################################
#Averaging between duplicate probes, pivoting longer, merging with patient ID, filtering for SHH patients only, filtering for 'intact' values, and pivoting wider before looping
Merge_DF <- readRDS("Finalized_Expression_Array.rds")
YBX_Panel <- Merge_DF[,-c(1:2)]

#Removing rows with no gene identifier - microarray controls and non-coding genes
YBX_Panel <- YBX_Panel[!(YBX_Panel$Gene_Name == '' | YBX_Panel$Gene_Name == '//' | YBX_Panel$Gene_Name == '//' | YBX_Panel$Gene_Name == 'NONCODE' | YBX_Panel$Gene_Name == '---' | YBX_Panel$Gene_Name == '///' | YBX_Panel$Gene_Name == 'ENSEMBL'),]

#Merging duplicate probes by mean
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
YBX_Cor_Results <- readRDS("C:/Users/12298/Desktop/Data_Analytics/Taylor_2017/YBX1_Cor_Results")
setwd("C:/Users/12298/Desktop/Data_Analytics/Taylor_2017")

cor(as.numeric(YBX_Cor_Test_2[17970,]), as.numeric(YBX_Cor_Test_2[16916,]))
plot(as.numeric(YBX_Cor_Test_2[17970,]), as.numeric(YBX_Cor_Test_2[16916,]), ylab="PARP1", xlab="YBX1", )

YBX1 <- YBX_Cor_Test_2[c(17970,16916),] %>% t() %>% as.data.frame()

ggplot(YBX1, aes(x=YBX1,y=PARP1)) + geom_point(pch=16, size=2, alpha=1, color = "#EF70DE") + ggtitle("YBX1-PARP1 Correlation in MB SHH Patients") + geom_smooth(method=lm, color="skyblue", fill="black", alpha=0.5) + geom_text(x=10.6, y=10.9, label="r2 = 0.6366", size=5) + theme_minimal(base_size = 16) + theme( legend.position = "none")

#End-----





