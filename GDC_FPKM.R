## ---------------------------
## Script name: GDC_FPKM.R
## Purpose of script: Do filtering and normalize readCount matrix in FPKM/UQ/TPM/CPM/RPKM
## Author: Dr. Nitish Kumar Mishra
## Date Created: 2021-03-02
## Copyright (c) Nitish Kumar Mishra, 2021
## Email: nitish.mishra@unmc.edu
## ---------------------------
setwd("C:/Users/nitis/Desktop/Davis RNAseq DEG")
library(dplyr)
gene_lengths <- read.csv("gencode.v37.table.txt", header = TRUE, sep = "\t")
rownames(gene_lengths) <- gene_lengths$Geneid
exp <- read.csv("STARreadContMatrix.txt", header = TRUE, sep="\t", row.names=1, 
                         stringsAsFactors = FALSE, check.names = FALSE)
colnames(exp) <- substr(colnames(exp), 1, nchar(colnames(exp))-1)# Remove last "_ " in colnames

########### Filtering Step/remove probes #############
Gene_filtering <- function(mat){
not_all_na <- function(x) any(!is.na(x))
mat <- mat %>%
  select(where(not_all_na))
mat <- subset(mat, apply(mat, 1, sum) != 0)
#mat <- subset(mat, apply(mat, 2, sum) != 0)# not_all_na already did this
mat <- mat[apply(mat,1,function(x) sum(x==0))< ncol(mat)*0.75,]##Remove genes which have over 25% zero
mat <- subset(mat, apply(mat, 1, sum) >= 10)
return(mat)}

########### Function for the FPKM & FPKM UQ ##########
GDC_FPKM <- function(mat, gene_lengths, method = "FPKM") 
{
  genes <- intersect(rownames(mat), gene_lengths$Geneid)
  gene_lengths <- gene_lengths[genes,]
  mat <- mat[genes,]
  protein_coding <- gene_lengths[grep("protein_coding", gene_lengths$Class),]
  RC_g = mat
  RC_pc <- colSums(mat[rownames(protein_coding),])
  #mat.pc <- mat[rownames(protein_coding),]
  eff_leng <- gene_lengths$Length
  names(eff_leng) <- gene_lengths$Geneid
  #eff_leng <- t(eff_leng)
  if (method == "FPKM") {
    fpkm <- do.call(cbind, lapply(1:ncol(RC_g), function(i) {
      (((RC_g[, i]) * 1e+09) / (eff_leng *RC_pc[i]))
    }))
    colnames(fpkm) <- colnames(RC_g)
    rownames(fpkm) <- rownames(RC_g)
    return(fpkm)}
  
  if(method=="FPKM-UQ") {
    uqs <- apply(RC_g, 2, quantile, 0.75)
    fpkm <- do.call(cbind, lapply(1:ncol(RC_g), function(i) {
      (((RC_g[, i]) * 1e+09) / (eff_leng * uqs[i]))
    }))
    colnames(fpkm) <- colnames(RC_g)
    rownames(fpkm) <- rownames(RC_g)
    return(fpkm)}
}
######################################################
####### This part is for getting common genes ########
common_Gene <- function(mat, geneLength){
  genes <- intersect(rownames(mat), geneLength$Geneid)
  gene_lengths <- geneLength[genes,]
  mat <- mat[genes,]
  return(list(Expression=mat, GeneLength=gene_lengths))}

######################################################
# https://gist.github.com/slowkow/6e34ccb4d1311b8fe62e
# https://www.reneshbedre.com/blog/expression_units.html
# This RPKM, TPM, CPM is same as Python 
rpkm <- function(counts, lengths) {
  rate <- (counts / lengths) * 1e9
  rate / sum(counts) }

tpm <- function(counts, lengths) {
  rate <- (counts / lengths)*1e3
  (rate / sum(rate, na.rm = TRUE)) * 1e6}

cpm <- function(counts, lengths)
{  (counts * 1e6) /sum(counts)}

######################################################
########### Filter and save normalized data ##########
ExpData <- Gene_filtering(mat = exp) ## Save after filtering
Common_List <- common_Gene(ExpData, gene_lengths)
rpkms <- apply(Common_List$Expression, 2, function(x) rpkm(x, Common_List$GeneLength$Length))
tpms <- apply(Common_List$Expression, 2, function(x) tpm(x, Common_List$GeneLength$Length))
cpms <- apply(Common_List$Expression, 2, function(x) cpm(x, Common_List$GeneLength$Length))
mat_FPKM <- as.data.frame(GDC_FPKM(mat = Common_List$Expression, gene_lengths = Common_List$GeneLength, method = "FPKM"))
mat_FPKM.UQ <- as.data.frame(GDC_FPKM(mat = Common_List$Expression, gene_lengths = Common_List$GeneLength, method = "FPKM-UQ"))

######################################################
### Mitochondrial (ChrM) have very high counts; exclude them for normalization
gene_lengths_ChrM <- gene_lengths[gene_lengths$Chromosome!="chrM",]
Common_List_ChrM <- common_Gene(ExpData, gene_lengths_ChrM)

rpkms.chrM <- apply(Common_List_ChrM$Expression, 2, function(x) rpkm(x, Common_List_ChrM$GeneLength$Length))
tpms.chrM <- apply(Common_List_ChrM$Expression, 2, function(x) tpm(x, Common_List_ChrM$GeneLength$Length))
cpms.chrM <- apply(Common_List_ChrM$Expression, 2, function(x) cpm(x, Common_List_ChrM$GeneLength$Length))

mat_FPKM.ChrM <- as.data.frame(GDC_FPKM(mat = Common_List_ChrM$Expression, gene_lengths = Common_List_ChrM$GeneLength, method = "FPKM"))
mat_FPKM.UQ.ChrM <- as.data.frame(GDC_FPKM(mat = Common_List_ChrM$Expression, gene_lengths = Common_List_ChrM$GeneLength, method = "FPKM-UQ"))

rm(common_Gene, Gene_filtering, cpm, GDC_FPKM, rpkm, tpm)
#################### Save the data ###################
save.image("NormalizedReads.RData")

#### Heatmap plot for each pathways ##################
######################################################
GENCODE.V37 <- gene_lengths_ChrM %>%
  select(Geneid, GeneSymbol, Class)%>%
  filter(stringr::str_detect(Class, " protein_coding"))

sampletable <- read.table("Complete_Samples_csv_file.txt", header=F, sep=",")
sampletable <- sampletable %>%
  select(Sample_Name=V1, Class=V4) %>%
  filter(!grepl("YOUNG_NonResp", Class)) %>%
  filter(!grepl("24_S9", Sample_Name)) %>%
  arrange(match(Class, c("YOUNG_Resp", "OLD_Highresp", "OLD_Nonresp"))) ## rearrange https://stackoverflow.com/questions/46129322/arranging-rows-in-custom-order-using-dplyr
  
rownames(sampletable) <- sampletable$Sample_Name

dataset1 <- rpkms.chrM[,sampletable$Sample_Name]
dataset1 <- as.data.frame(log2(dataset1+1))
dataset1$var <- matrixStats::rowVars(as.matrix(dataset1))
dataset1 <- merge(dataset1, GENCODE.V37, by=0)
dataset1$GeneSymbol <- gsub(" ", "", dataset1$GeneSymbol)
#rownames(dataset1) <- dataset1$GeneSymbol
pd <- sampletable$Class
ann_colors <- list(Group = c(YOUNG_Resp="black", OLD_Highresp = "green4", OLD_Nonresp="yellow3"))

##############################################
###### Change the name of pathway file #######
YAP_upregulated <- read.csv("Heatmap_for_R21/YAP_upregulated_genes.gmx", header = FALSE)[1]
Exp <- dataset1 %>% 
  filter(GeneSymbol %in% YAP_upregulated$V1)
rownames(Exp) <- Exp$GeneSymbol

library(tidyr)
Num_Gene <- 50
Exp <- Exp %>%
  slice_max(var, n=Num_Gene) %>%
  pivot_longer(grep("_S", names(Exp)),names_to = "Sample", values_to = "Expression") %>%
  mutate(Group=rep(pd, Num_Gene))

######## Change PDF file name for all pathways
tidyheatmap::tidy_heatmap(Exp, rows = GeneSymbol, columns = Sample,values = Expression, 
             scale = "row", annotation_col = c(Group), cluster_cols = FALSE,
             clustering_method = "ward.D2", clustering_distance_cols = "euclidean",
             #clustering_distance_cols = "manhattan",
             cluster_rows = TRUE, color_legend_n = 7,
             colors = c("blue4", "blue1","gray80","orange1","orange4"),
             annotation_colors = ann_colors, fontsize = 5, fontsize_row = 8,
             show_rownames = FALSE, fontsize_col = 5, angle_col = 90,
             height = 3, width = 3.5, filename = "YAP_upregulated.pdf"
             ) 


################# 
YAP_upregulated <- read.csv("Heatmap_for_R21/REACTOME_YAP1_AND_WWTR1_TAZ_STIMULATED_GENE_EXPRESSION .gmx", header = FALSE)[1]
Exp <- dataset1 %>% 
  filter(GeneSymbol %in% YAP_upregulated$V1)
rownames(Exp) <- Exp$GeneSymbol

library(tidyr)
Num_Gene <- 12
Exp <- Exp %>%
  slice_max(var, n=Num_Gene) %>%
  pivot_longer(grep("_S", names(Exp)),names_to = "Sample", values_to = "Expression") %>%
  mutate(Group=rep(pd, Num_Gene))

######## Change PDF file name for all pathways
tidyheatmap::tidy_heatmap(Exp, rows = GeneSymbol, columns = Sample,values = Expression, 
                          scale = "row", annotation_col = c(Group), cluster_cols = FALSE,
                          clustering_method = "ward.D2", clustering_distance_cols = "euclidean",
                          #clustering_distance_cols = "manhattan",
                          cluster_rows = TRUE, color_legend_n = 7,
                          colors = c("blue4", "blue1","gray80","orange1","orange4"),
                          annotation_colors = ann_colors, fontsize = 5, fontsize_row = 8,
                          show_rownames = FALSE, fontsize_col = 5, angle_col = 90,
                          height = 3, width = 3.5, filename = "REACTOME_YAP1_AND_WWTR1_TAZ.pdf"
) 
##################################################
save.image("Jitu Heatmap.RData")
