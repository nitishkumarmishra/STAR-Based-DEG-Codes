setwd("C:/Users/nitis/Desktop/HTseqReadCounts")
load("C:/Users/nitis/Desktop/HTseqReadCounts/HTseqRedCount.RData")
rm(list=setdiff(ls(), c("dataset", "PhenoData", "dataset1")))
################################################################
#https://jbengler.github.io/tidyheatmap/articles/tidyheatmap.html
library(tidyheatmap); library(dplyr); library(tidyverse); library(matrixStats)
############## GTF file ####################
gtf <- rtracklayer::import('gencode.v33.annotation.gtf')
gencode.v33.gtf <- as.data.frame(gtf)
gencode.v33.gtf.selected <- gencode.v33.gtf %>%
  filter(type=="gene") %>%
  rename(Chr=seqnames, Type=type,ENSG=gene_id, GeneType=gene_type, Symbol=gene_name) %>%
  select(Chr, Type, GeneType, ENSG, Symbol) 
############################################

diffExp <- read.csv("Group1 Vs Group2 and 3 DEG.csv", header = TRUE, sep = ",", row.names = 1, stringsAsFactors = FALSE, check.names = FALSE)
diffExp$variance <- rowVars(as.matrix(log2(diffExp[, grep("UR-", colnames(diffExp))]+1)), na.rm = TRUE)
diffExp <- diffExp[order(diffExp$pvalue, decreasing = FALSE),]
diffExp$ENSG <- diffExp$Gene

dataset1 <- diffExp %>%
  inner_join(gencode.v33.gtf.selected, by = "ENSG")


####################################
Num_Gene <- 100
Exp <- dataset1 %>%
  filter(GeneType %in% c("protein_coding", "lncRNA"))%>% 
  slice_min(pvalue, n=Num_Gene) %>%
  #slice(1:Num_Gene) %>%
  pivot_longer(cols = starts_with("UR"), names_to = "Sample", values_to = "Expression") %>%
  mutate(Group=rep(PhenoData$Group, Num_Gene)) %>%#rep(condition, 7); its seven gene so seventies of condition
  mutate(Expressio=log2(Expression + 1))

ann_colors <- list(Group = c("Group 1" = "red4", "Group 2" = "blue4", "Group 3"="green4"))

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
             #filename = "Heatmap1.pdf"
)
