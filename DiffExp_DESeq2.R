## On husker server
## Annotation:: gencode.v33.annotation.gtf [ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.chr_patch_hapl_scaff.annotation.gtf.gz]
## Reference:: GRCh38.p13.genome.fa [ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/GRCh38.p13.genome.fa.gz]
## STAR Version = 2.7.5a
## Step to prepapre readCounts are available at Evernote [Nagendra RNAseq]
setwd("C:/Users/nitis/Desktop/Davis RNAseq DEG")
library(DESeq2); library(dplyr)
count_matrix <- read.csv("STARreadContMatrix.txt", header = TRUE, sep="\t", row.names=1, 
                         stringsAsFactors = FALSE, check.names = FALSE)
colnames(count_matrix) <- substr(colnames(count_matrix), 1, nchar(colnames(count_matrix))-1)

sampletable <- read.table("Complete_Samples_csv_file.txt", header=F, sep=",")
sampletable <- sampletable %>%
  select(Sample_Name=V1, Class=V4) ## This will select and rename both at time
rownames(sampletable) <- sampletable$Sample_Name

dataset1 <- count_matrix[apply(count_matrix,1,function(x) sum(x==0))<ncol(count_matrix)*0.75,]##Remove genes which have over 25% zero
dataset1 <- dataset1[apply(dataset1,2,function(x) sum(x==0))<nrow(dataset1)*0.75,]##Remove samples which have over 25% zero
keep <- rowSums(dataset1) >= 10 ## Remove very low readcount genes.
dataset1 <- dataset1[keep,]

## OLD_Highresp Vs. OLD_Nonresp
PhenoData <- sampletable %>%
  filter(Class=="OLD_Highresp"| Class=="OLD_Nonresp")
#dataset1 <- dataset1[, rownames(PhenoData)]
#dataset_OLD_Resp_vs_OLD_Non <- dataset1[, c("45_S23", "4A_S5",  "13_S7",  "12_S6", "46_S24", "44_S22", "23_S8", "25_S10", "4_S4", "30_S12", "33_S13", "24_S9")]
countdata <- dataset1[,PhenoData$Sample_Name]
condition <- factor(ifelse(PhenoData$Class=="OLD_Highresp", "exp", "control"), levels = c("control", "exp"))
# Create a coldata frame and instantiate the DESeqDataSet. See ?DESeqDataSetFromMatrix
coldata <- data.frame(row.names=colnames(countdata), condition)
dds <- DESeqDataSetFromMatrix(countData=countdata, colData=coldata, design=~condition)
dds <- DESeq(dds); res <- results(dds)
table(res$padj<0.05)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=FALSE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "Gene"
#head(resdata)
resdata <- na.omit(resdata) ## Remove genes with Pvalue = NA
gencode.v37.gtf <- read.csv("gencode.v37.table.txt", header = TRUE, sep = "\t")
gencode.v37.gtf.selected <- gencode.v37.gtf %>%
  select(Gene=Geneid, Symbol=GeneSymbol, Chr=Chromosome, GeneType=Class)
resdata <- left_join(resdata,gencode.v37.gtf.selected, by ="Gene")
resdata <- resdata %>% 
  select(c("Gene","Symbol","Chr","GeneType"), everything()) ## Reoder result file
resdata_adjP_0.05 <- resdata %>%
  filter(stringr::str_detect(GeneType, "protein_coding") & padj <= 0.05)
 
##################################################
##################################################
save.image("DiffEx_DESeq2.RData")

