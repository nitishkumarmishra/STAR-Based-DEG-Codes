## ---------------------------
## Script name: GDC_FPKM.R
## Purpose of script: Do filtering and normalize readCount matrix in FPKM/UQ/TPM/CPM/RPKM
## Author: Dr. Nitish Kumar Mishra
## Date Created: 2021-03-02
## Copyright (c) Nitish Kumar Mishra, 2021
## Email: nitish.mishra@unmc.edu
## ---------------------------
setwd("C:/Users/nitis/Desktop/NAgendra STAR GeneCounts/ReadCount V37/")
library(dplyr)
library(limma)
## This part when rtracklayer was not working
#gene_lengths <- read.csv("gencode.v37.table.txt", header = TRUE, sep = "\t") 
gtf <- rtracklayer::import("gencode.v37.annotation.gtf")
gencode.v37.gtf <- as.data.frame(gtf)
gencode.v37.gtf.selected <- gencode.v37.gtf %>%
  filter(type=="gene") %>%
  rename(Chromosome= seqnames, Type= type, Geneid = gene_id, Class=gene_type, GeneSymbol= gene_name, Length=width, Start = start, End=end, Strand=strand) %>%
  select(Geneid, GeneSymbol, Chromosome, Start, End, Class,  Strand, Length) %>%
  filter(!grepl("chrM",Chromosome))
gene_lengths <- gencode.v37.gtf.selected  

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

#exp <- Gene_filtering(exp)
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

######################################################
############# Differential expression ################
sampletable <- as.data.frame(cbind(colnames(rpkms), 
                                   c(rep("DMSO", 3), rep("JOI", 3), rep("JOI_Plus_PAN", 3), rep("PAN", 3))))
sampletable <- sampletable %>%
  select(Sample_Name=V1, Class=V2) ## This will select and rename both at time
rownames(sampletable) <- sampletable$Sample_Name

ExpAnalysis <- function(dataset1=rpkms.chrM, case="PAN", control="DMSO", discard=NA ,pAdj=0.05, log2FC = 1)
{  
  PhenoData <- sampletable %>%
    filter(Class %in% case| Class==control)
  if(!is.na(discard))
  {
    PhenoData <- PhenoData[!PhenoData$Sample_Name%in%discard,]
    PhenoData <- PhenoData %>%
      filter(Class== case| Class==control)
  }
  
  countdata <- dataset1[,PhenoData$Sample_Name]
  condition <- factor(ifelse(PhenoData$Class %in% case, "exp", "control"), levels = c("control", "exp"))
  
  design <- model.matrix(~0+condition)
  colnames(design) <- c("control", "exp")

  cont_matrix <- makeContrasts(ControlVsExp = exp-control, levels=design)
  
  fit <- lmFit(countdata, design)
  fit_contrast <- contrasts.fit(fit, cont_matrix)
  fit_contrast <- eBayes(fit_contrast)
  top_genes <- topTable(fit_contrast, number = nrow(countdata), adjust = "BH", resort.by = "logFC",sort.by = "p", lfc = 1, p.value = 0.01)
  #top_genes <- topTableF(fit_contrast, number = nrow(countdata), adjust = "BH",sort.by = "F", lfc = 1, p.value = 0.01)
  resdata <- merge(as.data.frame(top_genes), as.data.frame(countdata), by="row.names", sort=FALSE)
  names(resdata)[1] <- "Geneid"
  resdata <- na.omit(resdata) ## Remove genes with Pvalue = NA
  resdata <- left_join(resdata,gencode.v37.gtf.selected, by ="Geneid")
  resdata <- resdata %>% 
    select(c("Geneid","GeneSymbol","Chromosome","Class"), everything()) ## Reoder result file
  
  resdata <- resdata %>%
    mutate(AvgCase=rowMeans(select(resdata, starts_with(case)))) %>%
    mutate(AvgControl=rowMeans(select(resdata, starts_with(control)))) %>%
    mutate(FoldChange=AvgCase/AvgControl)%>%
    mutate(log2FC_Actual=log2(FoldChange))
  resdata_adjP_0.05 <- resdata %>%
    filter(stringr::str_detect(Class, "protein_coding") & adj.P.Val <= pAdj & abs(logFC) >= log2FC)
  return(list(PhenoData=PhenoData, condition=condition, Result = resdata_adjP_0.05))
}

JOI_Vs_DMSO <- ExpAnalysis(dataset1 = cpms.chrM, case="JOI", control="DMSO", pAdj = 0.01, log2FC =  2)
PAN_Vs_DMSO <- ExpAnalysis(dataset1 = cpms.chrM, case="PAN", control="DMSO", pAdj = 0.01, log2FC = 2)
JOI_Vs_PAN <- ExpAnalysis(dataset1 = cpms.chrM, case="JOI", control="PAN", pAdj = 0.01, log2FC = 2)
JOI_PAN_Vs_DMSO <- ExpAnalysis(dataset1 = cpms.chrM, case=c("JOI","PAN"), control="DMSO", pAdj = 0.01, log2FC = 2)

JOI_Vs_DMSO <- ExpAnalysis(dataset1 = mat_FPKM.UQ.ChrM, case="JOI", control="DMSO", pAdj = 0.01, log2FC =  2)
PAN_Vs_DMSO <- ExpAnalysis(dataset1 = mat_FPKM.UQ.ChrM, case="PAN", control="DMSO", pAdj = 0.01, log2FC = 2)
JOI_Vs_PAN <- ExpAnalysis(dataset1 = mat_FPKM.UQ.ChrM, case="JOI", control="PAN", pAdj = 0.01, log2FC = 2)
JOI_PAN_Vs_DMSO <- ExpAnalysis(dataset1 = mat_FPKM.UQ.ChrM, case=c("JOI","PAN"), control="DMSO", pAdj = 0.01, log2FC = 2)

######################################################
######## Save analysis results in CSV files ##########
write.csv(JOI_Vs_DMSO$Result, "JOI_Vs_DMSO_TPM.csv")
write.csv(PAN_Vs_DMSO$Result, "PAN_Vs_DMSO_TPM.csv")
write.csv(JOI_Vs_PAN$Result, "JOI_Vs_PAN_TPM.csv")
write.csv(JOI_PAN_Vs_DMSO$Result,"JOI_PAN_Vs_DMSO_TPM.csv")

write.csv(JOI_Vs_DMSO$Result, "JOI_Vs_DMSO_CPM.csv")
write.csv(PAN_Vs_DMSO$Result, "PAN_Vs_DMSO_CPM.csv")
write.csv(JOI_Vs_PAN$Result, "JOI_Vs_PAN_CPM.csv")
write.csv(JOI_PAN_Vs_DMSO$Result,"JOI_PAN_Vs_DMSO_CPM.csv")


write.csv(JOI_Vs_DMSO$Result, "JOI_Vs_DMSO_RPKM.csv")
write.csv(PAN_Vs_DMSO$Result, "PAN_Vs_DMSO_RPKM.csv")
write.csv(JOI_Vs_PAN$Result, "JOI_Vs_PAN_RPKM.csv")
write.csv(JOI_PAN_Vs_DMSO$Result,"JOI_PAN_Vs_DMSO_RPKM.csv")


write.csv(JOI_Vs_DMSO$Result, "JOI_Vs_DMSO_FPKM.csv")
write.csv(PAN_Vs_DMSO$Result, "PAN_Vs_DMSO_FPKM.csv")
write.csv(JOI_Vs_PAN$Result, "JOI_Vs_PAN_FPKM.csv")
write.csv(JOI_PAN_Vs_DMSO$Result,"JOI_PAN_Vs_DMSO_FPKM.csv")


write.csv(JOI_Vs_DMSO$Result, "JOI_Vs_DMSO_FPKM-UQ.csv")
write.csv(PAN_Vs_DMSO$Result, "PAN_Vs_DMSO_FPKM-UQ.csv")
write.csv(JOI_Vs_PAN$Result, "JOI_Vs_PAN_FPKM-UQ.csv")
write.csv(JOI_PAN_Vs_DMSO$Result,"JOI_PAN_Vs_DMSO_FPKM-UQ.csv")

#################### Save the data ###################
save.image("NormalizedReads.RData")
######################################################