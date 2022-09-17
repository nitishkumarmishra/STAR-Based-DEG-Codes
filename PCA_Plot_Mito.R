##############################################
###### Change the name of pathway file #######
library(dplyr); library(plotly)
setwd("C:/Users/nitis/Desktop/Davis RNAseq DEG")
#load("NormalizedReads.RData")
load("Jitu Heatmap.RData")## Output workspace of GDC_FPKM.R

######################################################
GENCODE.V37 <- gene_lengths_ChrM %>%
  select(Geneid, GeneSymbol, Class)%>%
  filter(stringr::str_detect(Class, " protein_coding"))

sampletable <- read.table("Complete_Samples_csv_file.txt", header=F, sep=",")
sampletable <- sampletable %>%
  select(Sample_Name=V1, Class=V4) %>%
  #filter(!grepl("YOUNG_NonResp", Class)) %>%
  filter(!grepl("24_S9", Sample_Name)) %>%
  arrange(match(Class, c("YOUNG_Resp", "OLD_Highresp","YOUNG_NonResp", "OLD_Nonresp"))) ## rearrange https://stackoverflow.com/questions/46129322/arranging-rows-in-custom-order-using-dplyr

rownames(sampletable) <- sampletable$Sample_Name

dataset1 <- rpkms.chrM[,sampletable$Sample_Name]
dataset1 <- as.data.frame(log2(dataset1+1))
dataset1$var <- matrixStats::rowVars(as.matrix(dataset1))
dataset1 <- merge(dataset1, GENCODE.V37, by=0)
dataset1$GeneSymbol <- gsub(" ", "", dataset1$GeneSymbol)

######################################################

Mitichindial_Gene <- read.csv("Human.MitoCarta3.0 BroadInstitute.txt", header = TRUE, sep = "\t")
Exp <- dataset1 %>% 
  filter(GeneSymbol %in% Mitichindial_Gene$Symbol)
rownames(Exp) <- Exp$GeneSymbol

pd <- sampletable$Class
ann_colors <- list(Group = c(YOUNG_Resp="black", OLD_Highresp = "green4", YOUNG_NonResp="brown", OLD_Nonresp="yellow3"))

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
                          height = 3, width = 3.5, filename = "Mito top 50.pdf"
) 

#######################################################
#######################################################
library(plotly)

mat <- dataset1 %>% 
  filter(GeneSymbol %in% Mitichindial_Gene$Symbol)
rownames(mat) <- mat$GeneSymbol

#sampletable <- sampletable[colnames(mat),]
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

plotPCA3D(mat = mat, file = "3D-PCA-Plot_mito-500", ntop = 500)
#plotPCA3D(mat = mat, file = "3D-PCA-Plot_no-24_S9_rpkm", discard = TRUE, sample = "24_S9")

######### PCA plot autoplotly #########
library(autoplotly)
data <- t(mat)
data <- as.data.frame(merge(sampletable, data, by="row.names"))
autoplotly(prcomp(data[,4:ncol(data)]), data=data,
           colour = 'Class', frame = TRUE)
##########################################################
##########################################################
################## PCA plot ####################
library(data.table)
library(ggplot2)
Meth <- myNorm
#colnames(Meth) <- paste0(myLoad$pd$Slide, "_", myLoad$pd$Array)

Meth <- as.data.frame(mat)
pca <- prcomp(t(Meth), scale. = TRUE)
color <- as.factor(sampletable$Class)
PCi<-data.frame(pca$x,Sample=color)

ggplot(PCi,aes(x=PC1,y=PC2, col=Sample))+
  geom_point(size=4,alpha=0.9)+ #Size and alpha just for fun
  scale_color_manual(values = c("blue4","red4", "green4"))+ #your colors here
  theme_classic()
ggsave("2d-PCA plot mitochondial gene.pdf", dpi = 600, width = 5, height = 5)

################### 3D PCA plot ######################
groups <- levels(color)
colors <- c(rep("red4", 6), rep("blue4", 5), rep("green4", 6))
pca$pcolor <- colors

s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=1.2, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), angle = 55)
s3d.coords <- s3d$xyz.convert(pca$x[, 1]/100, pca$x[, 2]/100, pca$x[, 3]/100)
legend("top", inset = 0.05, bty = "n", cex = 1.2, groups, col = c("red4", "blue4", "green4"), pch = 19)
dev.print(pdf, '3d-PCA Plot mitochondial gene.pdf', width = 10, height = 10)


########### 3D-PCA plot with sample name #############
s3d <- scatterplot3d::scatterplot3d(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80, color = pca$pcolor, pch = 19, type = "h", lty.hplot = 0, cex.axis = 1.2, cex.lab = 1.2, cex.symbols=0.01, scale.y = 0.75, 
                                    xlab = paste("Comp 1: ", round(pca$sdev[1]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    ylab = paste("Comp 2: ", round(pca$sdev[2]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), 
                                    zlab = paste("Comp 3: ", round(pca$sdev[3]^2/sum(pca$sdev^2), 3)*100, "%", sep = ""), angle = 55)

s3d.coords <- s3d$xyz.convert(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80)
legend("top", inset = 0.05, bty = "n", cex = 1.2, groups, col = c("red4", "blue4", "green4", "black"), pch = 19)
#text(s3d$xyz.convert(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80),labels = substr(rownames(pca$x), 10, 17), cex= 0.6, col = pca$pcolor)
text(s3d$xyz.convert(pca$x[, 1]/80, pca$x[, 2]/80, pca$x[, 3]/80),labels = rownames(pca$x), cex= 0.6, col = pca$pcolor)
dev.print(pdf, '3d-PCA Plot With Name.pdf', width = 12, height = 12)


################### Z-score #########################
sampletable1 <- sampletable %>%
  filter(!grepl("YOUNG_Resp", Class)) %>%
  filter(!grepl("YOUNG_NonResp", Class)) %>%
  arrange(match(Class, c("OLD_Highresp", "OLD_Nonresp"))) ## rearrange https://stackoverflow.com/questions/46129322/arranging-rows-in-custom-order-using-dplyr

dataset2 <- dataset1
dataset2 = dataset2[!duplicated(dataset2$GeneSymbol),]
rownames(dataset2) <- dataset2$GeneSymbol
dataset2 <- dataset2[,sampletable1$Sample_Name]


sampletable2 <- sampletable1 %>%
  filter(!grepl("OLD_Nonresp", Class))

dataset2.OLD_Highresp <- dataset2[,sampletable2$Sample_Name]

sampletable2 <- sampletable1 %>%
  filter(!grepl("OLD_Highresp", Class))

dataset2.OLD_Nonresp <- dataset2[,sampletable2$Sample_Name]


Hippo_gene <- read.csv("Heatmap_for_R21/gene_sets.gmt", header = FALSE)

dataset2.OLD_Highresp.Hippo <- na.omit(dataset2.OLD_Highresp[Hippo_gene$V1,])
dataset2.OLD_Nonresp.Hippo <- na.omit(dataset2.OLD_Nonresp[Hippo_gene$V1,])

meanHippo_NonResp <- rowMeans(dataset2.OLD_Nonresp.Hippo)

#z_scores <- (dataset2.OLD_Highresp.Hippo - meanHippo_NonResp)/ matrixStats::rowSds(as.matrix(dataset2.OLD_Highresp.Hippo))

activation_status <- rowMeans(dataset2.OLD_Highresp.Hippo)/meanHippo_NonResp
sum(ifelse(activation_status > 1, + 1, -1))/sqrt(nrow(z_scores))


#######################################################
#######################################################
save.image("Davis_Mitochondrial_PCA.RData")

#######################################################
