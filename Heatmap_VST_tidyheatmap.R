## ---------------------------
## Script name: 
## Purpose of script:
## Author: Dr. Nitish Kumar Mishra
## Date Created: 2020-07-04
## Copyright (c) Nitish Kumar Mishra, 2020
## Email: nitish.mishra@unmc.edu
## ---------------------------
setwd("C:/Users/nitis/Desktop/HTseqReadCounts")

#################################################
pacman::p_unload(all)
library(dplyr)
############## GTF file ####################
gtf <- rtracklayer::import('gencode.v33.annotation.gtf')
gencode.v33.gtf <- as.data.frame(gtf)
gencode.v33.gtf.selected <- gencode.v33.gtf %>%
  filter(type=="gene") %>%
  rename(Chr=seqnames, Type=type, ENSG=gene_id, GeneType=gene_type, Symbol=gene_name) %>%
  select(Chr, Type, GeneType, ENSG, Symbol) 
#################################################
library(data.table)
Sample_list <-  list.files(path = ".",pattern=".tab")

for (file in Sample_list){
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    dataset <- read.table(file, header=TRUE, sep="\t",skip = 3, stringsAsFactors = FALSE)
    colnames(dataset) <- c("gene_id","Counts")
    dataset <- subset(dataset,select = c("gene_id","Counts"))
    setnames(dataset,old = c("Counts"),new =(file))  }
  # if the merged dataset does exist, append to it
  else {
    temp_dataset <- read.table(file, header=TRUE, sep="\t",skip = 3, stringsAsFactors = FALSE)
    colnames(temp_dataset) <- c("gene_id","Counts")
    temp_dataset <- subset(temp_dataset,select = c("gene_id","Counts"))
    setnames(temp_dataset,old = c("Counts"),new =(file))
    dataset<-merge(dataset, temp_dataset,by = 'gene_id')
    rm(temp_dataset)  }}

colnames(dataset) <- gsub("_R_ReadsPerGene.out.tab", "", colnames(dataset))
colnames(dataset) <- gsub("_-_|_", "-", colnames(dataset))

rownames(dataset) <- dataset$`gene-id`
dataset$`gene-id` <- NULL
colnames(dataset) <- substr(colnames(dataset), 1, 8)

##############################################

dataset1 <- dataset[, -grep("OBBM", colnames(dataset))]
dataset1 <- dataset1[apply(dataset1,1,function(x) sum(x==0))<ncol(dataset1)*0.75,]
dataset1 <- dataset1[apply(dataset1,2,function(x) sum(x==0))<nrow(dataset1)*0.75,]

PhenoData <- read.csv("PhenoData_NAS_Vs_Normal.txt", header = TRUE, sep = ",")



PhenoData$Group <- gsub(" ", "", PhenoData$Group)
condition <- factor(PhenoData$Group)
dataset1 <- dataset1[, as.character(PhenoData$Pheno)]

#PhenoData$Condition <- as.character(condition)

library(DESeq2)
countdata <- dataset1
# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
(coldata <- data.frame(row.names=colnames(countdata), condition))
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)


rld <- varianceStabilizingTransformation(dds)
Expression <- assay(rld)
variance <- matrixStats::rowVars(Expression)
Expression <- as.data.frame(Expression)
Expression$variance <- variance
Expression$ENSG <- rownames(Expression)

############################################

diffExp <- read.csv("Group1 Vs Group2 and 3 DEG.csv", header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
#diffExp$variance <- matrixStats::rowVars(as.matrix(log2(diffExp[, grep("UR-", colnames(diffExp))]+1)), na.rm = TRUE)
diffExp <- diffExp[order(diffExp$pvalue, decreasing = FALSE),]
diffExp$ENSG <- diffExp$Gene
diffExp <- diffExp %>%
  select(Gene, log2FoldChange, pvalue, padj, ENSG)

pacman::p_unload(all)
library(dplyr)
dataset1 <- Expression %>%
  inner_join(gencode.v33.gtf.selected, by = "ENSG")
dataset1 <- dataset1 %>%
  inner_join(diffExp, by = "ENSG") %>%
  rename(log2FC=log2FoldChange)%>%
  mutate(log2P= log2(pvalue), Direction=ifelse(log2FC > 0, "Up", "Down"))


####################################
library(tidyheatmap); library(tidyverse)
Num_Gene <- 100
Exp <- dataset1 %>%
  filter(GeneType %in% c("protein_coding", "lncRNA"))%>% 
  #filter(GeneType %in% c("protein_coding"))%>% 
  slice_min(pvalue, n=Num_Gene) %>%
  #slice_max(variance, n=Num_Gene) %>%
  #slice(1:Num_Gene) %>%
  pivot_longer(cols = starts_with("UR"), names_to = "Sample", values_to = "Expression") %>%
  mutate(Group=rep(PhenoData$Group, Num_Gene)) %>%#rep(condition, 7); its seven gene so seventies of condition
  mutate(Expressio=log2(Expression + 1))



ann_colors <- list(Group = c(Group1 = "red", Group2 = "blue", Group3="green"),
                   Direction=c(Up="red", Down="blue"),
                   log2FC= c("red","gray","blue"),
                   GeneType=c(protein_coding="red", lncRNA="blue"))


tidy_heatmap(Exp,
             rows = Symbol,
             columns = Sample,
             values = Expression,
             scale = "row",
             annotation_col = c(Group),
             annotation_row = c(GeneType, log2FC, Direction),
             cluster_cols = TRUE,
             clustering_method = "ward.D2",
             clustering_distance_cols = "manhattan",
             cluster_rows = TRUE,
             color_legend_n = 11,
             colors = c("green","gray","red"),
             annotation_colors = ann_colors,
             #angle_col = 315,
             height = 12,
             width = 14,
             filename = "Heatmap Group1 vs Group2 nad 3 Top 100.pdf"
)

########################################################
########################################################
PhenoData1 <- PhenoData[grep("Group3", PhenoData$Group, invert = TRUE),]
Expression1 <- Expression[,PhenoData1$Pheno]

variance <- matrixStats::rowVars(as.matrix(Expression1))
#Expression1 <- as.data.frame(Expression1)
Expression1$variance <- variance
Expression1$ENSG <- rownames(Expression)

diffExp <- read.csv("Group1 Vs Group2 DEG.csv", header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
#diffExp$variance <- matrixStats::rowVars(as.matrix(log2(diffExp[, grep("UR-", colnames(diffExp))]+1)), na.rm = TRUE)
diffExp <- diffExp[order(diffExp$pvalue, decreasing = FALSE),]
diffExp$ENSG <- diffExp$Gene
diffExp <- diffExp %>%
  select(Gene, log2FoldChange, pvalue, padj, ENSG)

pacman::p_unload(all)
library(dplyr)
dataset1 <- Expression1 %>%
  inner_join(gencode.v33.gtf.selected, by = "ENSG")
dataset1 <- dataset1 %>%
  inner_join(diffExp, by = "ENSG") %>%
  rename(log2FC=log2FoldChange)%>%
  mutate(log2P= log2(pvalue), Direction=ifelse(log2FC > 0, "Up", "Down"))

library(tidyheatmap); library(tidyverse)
Num_Gene <- 230
Exp <- dataset1 %>%
  filter(GeneType %in% c("protein_coding", "lncRNA"))%>% 
  #filter(GeneType %in% c("protein_coding"))%>% 
  slice_min(pvalue, n=Num_Gene) %>%
  #slice_max(variance, n=Num_Gene) %>%
  #slice(1:Num_Gene) %>%
  pivot_longer(cols = starts_with("UR"), names_to = "Sample", values_to = "Expression") %>%
  mutate(Group=rep(PhenoData1$Group, Num_Gene)) %>%#rep(condition, 7); its seven gene so seventies of condition
  mutate(Expressio=log2(Expression + 1))



ann_colors <- list(Group = c(Group1 = "red", Group2 = "blue"),
                   Direction=c(Up="red", Down="blue"),
                   log2FC= c("red","gray","blue"),
                   GeneType=c(protein_coding="red", lncRNA="blue"))


tidy_heatmap(Exp,
             rows = Symbol,
             columns = Sample,
             values = Expression,
             scale = "row",
             annotation_col = c(Group),
             annotation_row = c(GeneType, log2FC, Direction),
             cluster_cols = TRUE,
             clustering_method = "ward.D2",
             clustering_distance_cols = "manhattan",
             cluster_rows = TRUE,
             color_legend_n = 11,
             colors = c("green","gray","red"),
             annotation_colors = ann_colors,
             show_rownames = FALSE,
             fontsize = 10,
             #fontsize_row = 5,
             fontsize_col = 8,
             #angle_col = 315,
             height = 12,
             width = 14,
             filename = "Heatmap Group1 vs Group2 Top 230_1.pdf"
)
########################################################
# Group 1 and 2 Vs Group 3
########################################################
#PhenoData1 <- PhenoData[grep("Group3", PhenoData$Group, invert = TRUE),]
PhenoData1 <- PhenoData 
#PhenoData1$Group <- gsub("Group2", "Group1", PhenoData$Group)
Expression1 <- Expression[,PhenoData1$Pheno]
Expression1$ENSG <- rownames(Expression1)

diffExp <- read.csv("Group1 and 2 Vs Group3  DEG.csv", header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
#diffExp$variance <- matrixStats::rowVars(as.matrix(log2(diffExp[, grep("UR-", colnames(diffExp))]+1)), na.rm = TRUE)
diffExp <- diffExp[order(diffExp$pvalue, decreasing = FALSE),]
diffExp$ENSG <- diffExp$Gene
diffExp <- diffExp %>%
  select(Gene, log2FoldChange, pvalue, padj, ENSG)

pacman::p_unload(all)
library(dplyr)
dataset1 <- Expression1 %>%
  inner_join(gencode.v33.gtf.selected, by = "ENSG")
dataset1 <- dataset1 %>%
  inner_join(diffExp, by = "ENSG") %>%
  rename(log2FC=log2FoldChange)%>%
  mutate(log2P= log2(pvalue), Direction=ifelse(log2FC > 0, "Up", "Down"))

library(tidyheatmap); library(tidyverse)
Num_Gene <- 210
Exp <- dataset1 %>%
  filter(GeneType %in% c("protein_coding", "lncRNA"))%>% 
  #filter(GeneType %in% c("protein_coding"))%>% 
  slice_min(pvalue, n=Num_Gene) %>%
  #slice_max(variance, n=Num_Gene) %>%
  #slice(1:Num_Gene) %>%
  pivot_longer(cols = starts_with("UR"), names_to = "Sample", values_to = "Expression") %>%
  mutate(Group=rep(PhenoData1$Group, Num_Gene)) %>%#rep(condition, 7); its seven gene so seventies of condition
  mutate(Expressio=log2(Expression + 1))



ann_colors <- list(Group = c(Group1 = "red", Group2 = "blue", Group3="green"),
                   Direction=c(Up="red", Down="blue"),
                   log2FC= c("red","gray","blue"),
                   GeneType=c(protein_coding="red", lncRNA="blue"))


tidy_heatmap(Exp,
             rows = Symbol,
             columns = Sample,
             values = Expression,
             scale = "row",
             annotation_col = c(Group),
             annotation_row = c(GeneType, log2FC, Direction),
             cluster_cols = TRUE,
             clustering_method = "ward.D2",
             clustering_distance_cols = "manhattan",
             cluster_rows = TRUE,
             color_legend_n = 11,
             colors = c("green","gray","red"),
             annotation_colors = ann_colors,
             show_rownames = FALSE,
             fontsize = 10,
             #fontsize_row = 5,
             fontsize_col = 8,
             #angle_col = 315,
             height = 8,
             width = 10,
             filename = "Heatmap Group1 and 2 vs Group3 Top 210.pdf"
)

########################################################
# Group 1 Vs Group 3
########################################################
#PhenoData1 <- PhenoData[grep("Group3", PhenoData$Group, invert = TRUE),]
PhenoData1 <- PhenoData 
PhenoData1 <- PhenoData1[grep("Group2", PhenoData1$Group, invert = TRUE),]
Expression1 <- Expression[,PhenoData1$Pheno]
Expression1$ENSG <- rownames(Expression1)

diffExp <- read.csv("Group1 Vs Group3 DEG.csv", header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
#diffExp$variance <- matrixStats::rowVars(as.matrix(log2(diffExp[, grep("UR-", colnames(diffExp))]+1)), na.rm = TRUE)
diffExp <- diffExp[order(diffExp$pvalue, decreasing = FALSE),]
diffExp$ENSG <- diffExp$Gene
diffExp <- diffExp %>%
  select(Gene, log2FoldChange, pvalue, padj, ENSG)

pacman::p_unload(all)
library(dplyr)
dataset1 <- Expression1 %>%
  inner_join(gencode.v33.gtf.selected, by = "ENSG")
dataset1 <- dataset1 %>%
  inner_join(diffExp, by = "ENSG") %>%
  rename(log2FC=log2FoldChange)%>%
  mutate(log2P= log2(pvalue), Direction=ifelse(log2FC > 0, "Up", "Down"))

library(tidyheatmap); library(tidyverse)
Num_Gene <- 251
Exp <- dataset1 %>%
  filter(GeneType %in% c("protein_coding", "lncRNA"))%>% 
  #filter(GeneType %in% c("protein_coding"))%>% 
  slice_min(pvalue, n=Num_Gene) %>%
  #slice_max(variance, n=Num_Gene) %>%
  #slice(1:Num_Gene) %>%
  pivot_longer(cols = starts_with("UR"), names_to = "Sample", values_to = "Expression") %>%
  mutate(Group=rep(PhenoData1$Group, Num_Gene)) %>%#rep(condition, 7); its seven gene so seventies of condition
  mutate(Expressio=log2(Expression + 1))



ann_colors <- list(Group = c(Group1 = "red", Group3="green"),
                   Direction=c(Up="red", Down="blue"),
                   log2FC= c("red","gray","blue"),
                   GeneType=c(protein_coding="red", lncRNA="blue"))


tidy_heatmap(Exp,
             rows = Symbol,
             columns = Sample,
             values = Expression,
             scale = "row",
             annotation_col = c(Group),
             annotation_row = c(GeneType, log2FC, Direction),
             cluster_cols = TRUE,
             clustering_method = "ward.D2",
             clustering_distance_cols = "manhattan",
             cluster_rows = TRUE,
             color_legend_n = 11,
             colors = c("green","gray","red"),
             annotation_colors = ann_colors,
             show_rownames = FALSE,
             fontsize = 10,
             fontsize_row = 5,
             fontsize_col = 6,
             #angle_col = 315,
             height = 8,
             width = 10,
             filename = "Heatmap Group1 vs Group3 Top 215.pdf"
)

save.image("Heatmap_VST_tidyheatmap.RData")
