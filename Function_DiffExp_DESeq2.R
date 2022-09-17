## ---------------------------
## Script name: Function_DiffExp_DESeq2.R
## Purpose of script: Run DESeq2 based DEG analysis function by using function "ExpAnalysis"
## Author: Dr. Nitish Kumar Mishra
## Date Created: 2021-02-18
## Copyright (c) Nitish Kumar Mishra, 2021
## Email: nitish.mishra@unmc.edu
## ---------------------------
## On husker server
## Annotation:: gencode.v33.annotation.gtf [ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.chr_patch_hapl_scaff.annotation.gtf.gz]
## Reference:: GRCh38.p13.genome.fa [ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.p13.genome.fa.gz]
## STAR Version = 2.7.5a
## Step to prepapre readCounts are available at Evernote [Nagendra RNAseq]
############ Command to run the code ############
# OLD_Highresp_Vs_YOUNG_NonResp <- ExpAnalysis(case="OLD_Highresp", control="YOUNG_NonResp", pAdj = 0.05)
setwd("C:/Users/nitis/Desktop/Davis RNAseq DEG")
library(DESeq2); library(dplyr)
count_matrix <- read.csv("STARreadContMatrix.txt", header = TRUE, sep="\t", row.names=1, 
                         stringsAsFactors = FALSE, check.names = FALSE)
colnames(count_matrix) <- substr(colnames(count_matrix), 1, nchar(colnames(count_matrix))-1)
sampletable <- read.table("Complete_Samples_csv_file.txt", header=F, sep=",")
sampletable <- sampletable %>%
  select(Sample_Name=V1, Class=V4) ## This will select and rename both at time
rownames(sampletable) <- sampletable$Sample_Name
#### Filter out genes with very low expression
dataset1 <- count_matrix[apply(count_matrix,1,function(x) sum(x==0))<ncol(count_matrix)*0.75,]##Remove genes which have over 25% zero
dataset1 <- dataset1[apply(dataset1,2,function(x) sum(x==0))<nrow(dataset1)*0.75,]##Remove samples which have over 25% zero
keep <- rowSums(dataset1) >= 10 ## Remove very low readcount genes.
dataset1 <- dataset1[keep,]
gencode.v37.gtf <- read.csv("gencode.v37.table.txt", header = TRUE, sep = "\t")
gencode.v37.gtf.selected <- gencode.v37.gtf %>%
  select(Gene=Geneid, Symbol=GeneSymbol, Chr=Chromosome, GeneType=Class)

gencode.v37.gtf.selected <- gencode.v37.gtf.selected[gencode.v37.gtf.selected$Chr!="chrM",]

ExpAnalysis <- function(case="OLD_Highresp", control="OLD_Nonresp", discard=NA ,pAdj=0.05)
{  
  
  PhenoData <- sampletable %>%
    filter(Class== case| Class==control)
  if(!is.na(discard))
  {
    PhenoData <- PhenoData[!PhenoData$Sample_Name%in%discard,]
    PhenoData <- PhenoData %>%
      filter(Class== case| Class==control)
  }
  
  countdata <- dataset1[,PhenoData$Sample_Name]
  condition <- factor(ifelse(PhenoData$Class==case, "exp", "control"), levels = c("control", "exp"))
  coldata <- data.frame(row.names=colnames(countdata), condition)
  dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
  dds <- DESeq(dds);  res <- results(dds)
  resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=FALSE)), by="row.names", sort=FALSE)
  names(resdata)[1] <- "Gene"
  resdata <- na.omit(resdata) ## Remove genes with Pvalue = NA
  resdata <- left_join(resdata,gencode.v37.gtf.selected, by ="Gene")
  resdata <- resdata %>% 
    select(c("Gene","Symbol","Chr","GeneType"), everything()) ## Reoder result file
  resdata_adjP_0.05 <- resdata %>%
    filter(stringr::str_detect(GeneType, "protein_coding") & padj <= pAdj)
  return(list(PhenoData=PhenoData, condition=condition, Result = resdata_adjP_0.05))
}

OLD_Highresp_Vs_OLD_Nonresp <- ExpAnalysis(case="OLD_Highresp", control="OLD_Nonresp", pAdj = 0.05)
OLD_Highresp_Vs_YOUNG_NonResp <- ExpAnalysis(case="OLD_Highresp", control="YOUNG_NonResp", pAdj = 0.05)
OLD_Highresp_Vs_YOUNG_Resp <- ExpAnalysis(case="OLD_Highresp", control="YOUNG_Resp", pAdj = 0.05)
OLD_Nonresp_Vs_YOUNG_NonResp <- ExpAnalysis(case="OLD_Nonresp", control="YOUNG_NonResp", pAdj = 0.05)
OLD_Nonresp_Vs_YOUNG_Resp <- ExpAnalysis(case="OLD_Nonresp", control="YOUNG_Resp", pAdj = 0.05)
YOUNG_Resp_Vs_YOUNG_NonResp <- ExpAnalysis(case="YOUNG_Resp", control="YOUNG_NonResp", pAdj = 0.05)


OLD_Highresp_Vs_OLD_Nonresp <- ExpAnalysis(case="OLD_Highresp", control="OLD_Nonresp", pAdj = NA, discard ="24_S9")


keyvals <- ifelse(
  res$log2FoldChange < -1.5, 'royalblue',
  ifelse(res$log2FoldChange > 1.5, 'gold',
         'black'))
keyvals[is.na(keyvals)] <- 'black'
names(keyvals)[keyvals == 'gold'] <- 'high'
names(keyvals)[keyvals == 'black'] <- 'mid'
names(keyvals)[keyvals == 'royalblue'] <- 'low'

EnhancedVolcano(OLD_Highresp_Vs_OLD_Nonresp$Result, lab = OLD_Highresp_Vs_OLD_Nonresp$Result$Symbol, 
                x = 'log2FoldChange', y = 'padj', ylab = bquote(~-Log[10] ~ italic(adjP)),
                xlim = c(-5.5, 5.5),
                pCutoff = 0.05, FCcutoff = 1.5, 
                title = 'Old responder Vs Old non-responder',
                subtitle = 'Differential expression',
                col=c('black', 'black', 'black', 'red'),
                pointSize = c(ifelse(abs(OLD_Highresp_Vs_OLD_Nonresp$Result$log2FoldChange)>2 & OLD_Highresp_Vs_OLD_Nonresp$Result$padj < 0.0001, 3, 1)))
##############################################################
save.image("Function_DESeq2.RData")
##############################################################
## Sample to remove 45_S23, 41_S21, 39_S19, 37_S17

OLD_Highresp_Vs_OLD_Nonresp <- ExpAnalysis(case="OLD_Highresp", control="OLD_Nonresp", pAdj = 0.05, discard = c("45_S23", "41_S21", "39_S19", "37_S17"))
OLD_Highresp_Vs_YOUNG_NonResp <- ExpAnalysis(case="OLD_Highresp", control="YOUNG_NonResp", pAdj = 0.05, discard = c("45_S23", "41_S21", "39_S19", "37_S17"))
OLD_Highresp_Vs_YOUNG_Resp <- ExpAnalysis(case="OLD_Highresp", control="YOUNG_Resp", pAdj = 0.05, discard = c("45_S23", "41_S21", "39_S19", "37_S17"))
OLD_Nonresp_Vs_YOUNG_NonResp <- ExpAnalysis(case="OLD_Nonresp", control="YOUNG_NonResp", pAdj = 0.05, discard = c("45_S23", "41_S21", "39_S19", "37_S17"))
OLD_Nonresp_Vs_YOUNG_Resp <- ExpAnalysis(case="OLD_Nonresp", control="YOUNG_Resp", pAdj = 0.05, discard = c("45_S23", "41_S21", "39_S19", "37_S17"))
YOUNG_Resp_Vs_YOUNG_NonResp <- ExpAnalysis(case="YOUNG_Resp", control="YOUNG_NonResp", pAdj = 0.05, discard = c("45_S23", "41_S21", "39_S19", "37_S17"))

save.image("DESeq2_Discard.RData")
#############################################
########## PCA plot from autoplotly #########
library(autoplotly)
data <- t(count_matrix)
data <- as.data.frame(merge(sampletable, data, by="row.names"))
autoplotly(prcomp(data[,4:29027]), data=data,
           colour = 'Class', frame = TRUE)

