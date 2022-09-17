### Differential Gene Expression analysis DESeq2
setwd("C:/Users/nitis/Desktop/NAgendra STAR GeneCounts/ReadCount V37/")
library(DESeq2); library(dplyr)
count_matrix <- read.csv("STARreadContMatrix.txt", header = TRUE, sep="\t", row.names=1, 
                         stringsAsFactors = FALSE, check.names = FALSE)
colnames(count_matrix) <- substr(colnames(count_matrix), 1, nchar(colnames(count_matrix))-1)
sampletable <- as.data.frame(cbind(colnames(count_matrix), 
                                   c(rep("DMSO", 3), rep("JOI", 3), rep("JOI_PAN", 3), rep("PAN", 3))))
sampletable <- sampletable %>%
  select(Sample_Name=V1, Class=V2) ## This will select and rename both at time
rownames(sampletable) <- sampletable$Sample_Name
#### Filter out genes with very low expression
dataset1 <- count_matrix[apply(count_matrix,1,function(x) sum(x==0))<ncol(count_matrix)*0.75,]##Remove genes which have over 25% zero
dataset1 <- dataset1[apply(dataset1,2,function(x) sum(x==0))<nrow(dataset1)*0.75,]##Remove samples which have over 25% zero
keep <- rowSums(dataset1) >= 10 ## Remove very low readcount genes.
dataset1 <- dataset1[keep,]
# gencode.v37.gtf <- read.csv("gencode.v37.table.txt", header = TRUE, sep = "\t")
# gencode.v37.gtf.selected <- gencode.v37.gtf %>%
#   select(Gene=Geneid, Symbol=GeneSymbol, Chr=Chromosome, GeneType=Class)%>%
#   filter(!grepl("chrM",Chr)) ## Remove gene mapped on chrM, redas are unequally mapper on mitochondia.

gtf <- rtracklayer::import("gencode.v37.annotation.gtf")
gencode.v37.gtf <- as.data.frame(gtf)
gencode.v37.gtf.selected <- gencode.v37.gtf %>%
  filter(type=="gene") %>%
  rename(Chr= seqnames, Type= type, Gene = gene_id, GeneType=gene_type, Symbol= gene_name, Length=width, Start = start, End=end, Strand=strand) %>%
  select(Gene, Symbol, Chr, Start, End, GeneType,  Strand, Length) %>%
  filter(!grepl("chrM",Chr))

ExpAnalysis <- function(case="PAN", control="DMSO", discard=NA ,pAdj=0.05, log2FC = 1)
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
    filter(stringr::str_detect(GeneType, "protein_coding") & padj <= pAdj & abs(log2FoldChange) >= log2FC)
  return(list(PhenoData=PhenoData, condition=condition, Result = resdata_adjP_0.05))
}

JOI_Vs_DMSO <- ExpAnalysis(case="JOI", control="DMSO", pAdj = 0.01, log2FC =  1)
PAN_Vs_DMSO <- ExpAnalysis(case="PAN", control="DMSO", pAdj = 0.01, log2FC = 1)
JOI_Vs_PAN <- ExpAnalysis(case="JOI", control="PAN", pAdj = 0.01, log2FC = 1)
JOI_PAN_Vs_DMSO <- ExpAnalysis(case=c("JOI","PAN"), control="DMSO", pAdj = 0.01, log2FC = 1)

###############################################
MYC <- count_matrix[rownames(count_matrix)=="ENSG00000136997.21",] ### ENSG00000136997 is MYC
MYC <- log2(MYC)
t.test(MYC[1:3], MYC[4:6])

##############################################################
write.csv(JOI_Vs_DMSO$Result, "JOI_Vs_DMSO_DESeq2.csv")
write.csv(PAN_Vs_DMSO$Result, "PAN_Vs_DMSO_DESeq2.csv")
write.csv(JOI_Vs_PAN$Result, "JOI_Vs_PAN_DESeq2.csv")
write.csv(JOI_PAN_Vs_DMSO$Result,"JOI_PAN_Vs_DMSO_DESeq2.csv")

##############################################################
save.image("DESeq2.RData")

