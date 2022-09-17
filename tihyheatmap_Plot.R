setwd("C:/Users/nitis/Desktop/HTseqReadCounts")
load("C:/Users/nitis/Desktop/HTseqReadCounts/HTseqRedCount.RData")
rm(list=setdiff(ls(), c("dataset", "PhenoData")))
################################################################
#https://jbengler.github.io/tidyheatmap/articles/tidyheatmap.html
library(tidyheatmap)
library(dplyr)
library(tidyverse)
library(matrixStats)
## DEG matrix is below, Excel sheet from Radha. Usea datapasta to make DEG as data.frame.
## 1.RNA seq-data_Group 1 vs Group 3-FC 1.5 with gene names to Nitish.xlsx
DEG <- data.frame(stringsAsFactors = FALSE, check.names = FALSE,
  ENSG = c("ENSG00000140465","ENSG00000084674","ENSG00000089169","ENSG00000179915", "ENSG00000242512","ENSG00000280660","ENSG00000144406"),
  GENE = c("CYP1A1","APOB","RPH3A","NRXN1","LINC01206","AL157396.1","UNC80"),
  baseMean = c(22.54392485,6.51709883,2.142333028,3.54332457,2.203491505,1.684918286,3.807731829),
  log2Fold.Change = c(4.299510273,-1.553326618,-2.847450134,-2.928603958,-2.979048942,-2.586232609, -2.573489481),
  lfcSE = c(0.66584089,0.292901663, 0.544347155,0.572837437,0.608721899,0.531798081,0.539380015),
  stat = c(6.457263794,-5.303235912,-5.230945191,-5.1124521,-4.89394081,-4.863185301,-4.771199172),
  `p-value` = c(1.07e-10,1.14e-07,1.69e-07, 3.18e-07,9.88e-07,1.16e-06,1.83e-06),
  `p-adj` = c(3.84e-06,0.002027288, 0.002027288,0.002867041,0.006225695,0.006225695,0.007338113))
#################################################################
PhenoData <- PhenoData %>%
  filter(Group=="Group 1"| Group=="Group 2") 

dataset1 <- dataset[, -grep("OBBM", colnames(dataset))]
rownames(dataset1) <- dataset1$`gene-id`; dataset1$`gene-id` <- NULL
dataset1 <- dataset1[apply(dataset1,1,function(x) sum(x==0))<ncol(dataset1)*0.75,]##Remove genes which have over 25% zero
dataset1 <- dataset1[apply(dataset1,2,function(x) sum(x==0))<nrow(dataset1)*0.75,]##Remove samples which have over 25% zero
keep <- rowSums(dataset1) >= 10 ## Remove very low readcount genes.
dataset1 <- dataset1[keep,]
colnames(dataset1) <- gsub("-Group-[0-9]", "", colnames(dataset1))
dataset1 <- dataset1[,PhenoData$Pheno]
dataset1 <- log2(dataset1+1)
#dataset1$variance <- rowVars(as.matrix(log2(dataset1+1)), na.rm = TRUE)
dataset1$variance <- rowVars(as.matrix(dataset1), na.rm = TRUE)
dataset1 <- dataset1[order(dataset1$variance, decreasing = TRUE),]
dataset1$ENSG = rownames(dataset1)

############## GTF file ####################
gtf <- rtracklayer::import('gencode.v33.annotation.gtf')
gencode.v33.gtf <- as.data.frame(gtf)
gencode.v33.gtf.selected <- gencode.v33.gtf %>%
  filter(type=="gene") %>%
  rename(Chr=seqnames, Type=type,ENSG=gene_id, GeneType=gene_type, Symbol=gene_name) %>%
  select(Chr, Type, GeneType, ENSG, Symbol) 
############################################
dataset1 <- dataset1 %>%
  inner_join(gencode.v33.gtf.selected, by = "ENSG")


#PhenoData <- PhenoData %>%
#  filter(Group=="Group 1"| Group=="Group 2") 
#dataset1 <- dataset1[, as.character(PhenoData$Pheno)]
#dataset1$ENSG <- sapply(strsplit(dataset1$ENSG, split = "[.]"),'[', 1)

#rowVars(as.matrix(dataset1), na.rm = TRUE)
#sapply(strsplit(head(rownames(dataset1)), split = "[.]"),'[', 1)
###############################################################
#Exp <- dataset1
#Exp <- dataset1[DEG$ENSG,]
#rownames(Exp) <- sapply(strsplit(rownames(Exp), split = "[.]"),'[', 1)
#Exp <- as.data.frame(t(Exp))
#Exp$Group <- c(rep("Group1", 32), rep("Group2", 32))
#Exp <- as.data.frame(t(Exp))
# rownames(DEG) <- DEG$GENE
# rownames(Exp) <- DEG$GENE
# DEG <- DEG %>%
#   mutate(Direction = ifelse(log2Fold.Change >0, "Up", "Down")) %>%
#   select(GENE,"p-adj", "Direction")
# Exp <-   merge(Exp, DEG, by.x = "row.names", by.y = "GENE")
# 
# Exp <- Exp %>%
#   rename(ENSG="Row.names") 
#   
# Exp <- Exp %>%
#   pivot_longer(cols = starts_with("UR"), names_to = "Sample", values_to = "Expression") %>%
#   mutate(Group=rep(c(rep("Group1", 32), rep("Group2", 32)),7))#rep(condition, 7); its seven gene so seventies of condition
# ############################################################
Num_Gene <- 500
Exp <- dataset1 %>%
  filter(GeneType %in% c("protein_coding", "lncRNA"))%>% 
  slice(1:Num_Gene) %>%
  pivot_longer(cols = starts_with("UR"), names_to = "Sample", values_to = "Expression") %>%
  mutate(Group=rep(PhenoData$Group, Num_Gene))#rep(condition, 7); its seven gene so seventies of condition


ann_colors <- list(    Group = c("Group 1" = "red4", "Group 2" = "blue4"))

tidy_heatmap(Exp,
             rows = Symbol,
             columns = Sample,
             values = Expression,
             scale = "row",
             annotation_col = c(Group),
             annotation_row = c(GeneType, Chr),
             cluster_cols = TRUE,
             clustering_method = "ward.D2",
             clustering_distance_cols = "manhattan",
             cluster_rows = TRUE,
             color_legend_n = 7,
             colors = c("red4","gray","green4"),
             annotation_colors = ann_colors,
             #angle_col = 315,
             height = 12,
             width = 14,
             #filename = "Heatmap.pdf"
             )
## Tidy z-score normalization
iris %>%
  effectsize::standardize() %>% 
  head()

