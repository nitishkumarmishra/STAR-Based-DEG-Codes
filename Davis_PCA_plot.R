#####
setwd("C:/Users/nitis/Desktop/Davis RNAseq DEG")
library(dplyr); library(plotly)
load("NormalizedReads.RData")

sampletable <- read.table("Complete_Samples_csv_file.txt", header=F, sep=",")
sampletable <- sampletable %>%
  select(Sample_Name=V1, Class=V4) %>%
  filter(!grepl("YOUNG_NonResp", Class))

rownames(sampletable) <- sampletable$Sample_Name

#sampletable <- sampletable[colnames(mat),]
mat <- as.data.frame(rpkms.chrM)
mat <- mat[,sampletable$Sample_Name]
### 3D PCA PLOT
## discard=TRUE, then in sample give the list of samples
plotPCA3D <- function (mat, filter=T, discard=TRUE, sample=NA, GeneSelect=TRUE, ntop = 500, file="3D-PCA-Plot"){
  #mat <- as.matrix(mat)
  not_all_na <- function(x) any(!is.na(x))
  mat <- mat %>%
    select(where(not_all_na))
  mat <- subset(mat, apply(mat, 1, sum) != 0)
  #mat <- subset(mat, apply(mat, 2, sum) != 0)# not_all_na already did this
  if(filter==T){
    mat <- mat[apply(mat,1,function(x) sum(x==0))< ncol(mat)*0.75,]##Remove genes which have over 25% zero
    mat <- subset(mat, apply(mat, 1, sum) >= 10)  }
  
  if((discard==TRUE) && (!is.na(sample)))
  {
    mat <- mat[, !names(mat)%in%sample]
    #sampletable <- sampletable[colnames(mat),]
    sampletable <- sampletable[!rownames(sampletable)%in%sample,]
  }
  rv <- matrixStats::rowVars(as.matrix(mat))
  if(GeneSelect==TRUE && (!is.na(ntop)))
  { 
    select <- order(rv, decreasing = TRUE)
    select <- select[1:ntop]
    mat <- mat[select,] 
  }
  #select <- order(rv, decreasing = TRUE)[seq_len(min(ntop, length(rv)))]
  pca <- prcomp(t(mat), scale. = TRUE)
  percentVar <- round(pca$sdev^2/sum(pca$sdev^2),3)*100
  group <- sampletable$Class
  d <- data.frame(PC1 = pca$x[, 1], PC2 = pca$x[, 2], PC3 = pca$x[, 3],
                  group = group, name = colnames(mat))
  message("Generating plotly plot")
  fig <- plotly::plot_ly(data = d, x = ~PC1, y = ~PC2, z = ~PC3, color = group, text = rownames(d))
  fig <- fig %>% add_markers()
  fig <- fig %>% layout(
    title = "Layout options in a 3d scatter plot",
    scene = list(
      xaxis = list(title = paste0("Comp 1: ", percentVar[1], "%", sep = "")),
      yaxis = list(title =  paste0("Comp 2: ", percentVar[2], "%", sep = "")),
      zaxis = list(title = paste0("Comp 3: ", percentVar[3], "%", sep = ""))
    ))
  #return(fig)
  htmlwidgets::saveWidget(as_widget(fig), paste0(file, ".html"))
}

plotPCA3D(mat = mat, file = "3D-PCA-Plot_rpkm")
plotPCA3D(mat = mat, file = "3D-PCA-Plot_no-24_S9_rpkm", discard = TRUE, sample = "24_S9")

######### PCA plot autoplotly #########
library(autoplotly)
data <- t(mat)
data <- as.data.frame(merge(sampletable, data, by="row.names"))
autoplotly(prcomp(data[,4:28159]), data=data,
           colour = 'Class', frame = TRUE)

#########################
save("PCA_Plot.RData")