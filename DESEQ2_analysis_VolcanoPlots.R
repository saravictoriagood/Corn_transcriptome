####################################################################################################
############################## HEATMAP ##################################################
####################################################################################################

## heatmap
# Extract gene expression values using significant genes
sig_gene_expression=gene_expression[rownames(gene_expression) %in% sig_tn_de$id,]
#remove control + tube columns
sig_gene_expression=sig_gene_expression[,-c(7:8)]


head(gene_expression)
head(sig_gene_expression)


# for pheatmap function, column names and row names of data and pdata mush be identical; change the row names
rownames(phenodata)=phenodata[,1]

# remove the id column
phenotype_table=subset(phenodata, select = -c(id) )

# change the colnames to match with the sample names
colnames(sig_gene_expression)=row.names(phenotype_table)

library(pheatmap)
pheatmap(as.matrix(sig_gene_expression), 
         scale = "row", 
         clustering_distance_rows = "correlation", 
         clustering_method = "complete",
         annotation_col = phenotype_table , 
         main="Significant genes",
         fontsize_col=14, 
         fontsize_row = 6,
         color = c("green","red"),
)


###############################################################################################################
##################################   PCA PLOTS ####################################################
###############################################################################################################
dev.off()
# PCA = Principal Component Analysis (PCA)
# A PCA plot shows clusters of samples based on their similarity.

## Draw PCA plot
# transpose the data and compute principal components
pca_data=prcomp(t(sig_gene_expression))


# Calculate PCA component percentages
pca_data_perc=round(100*pca_data$sdev^2/sum(pca_data$sdev^2),1)


# Extract 1 and 2 principle components and create a data frame 
# with sample names, first and second principal components and group information
df_pca_data = data.frame(PC1 = pca_data$x[,1], 
                         PC2 = pca_data$x[,2], 
                         sample = colnames(sig_gene_expression),
                         #condition = c("cont2","1451_1","1451_2","cont1","1451_3","cont3")
                         condition = c("control","1451","1451","control","1451","control")
)
)

head(sig_gene_expression)


## We will use ggplot2 and ggrepel (Repel overlapping text labels away from each other) packages to draw the PCA plots.

install.packages("ggrepel")

library(ggplot2)
library(ggrepel)

# color by sample
ggplot(df_pca_data, 
       aes(PC1,PC2, 
           color = sample))+
  geom_point(size=8)+
  labs(x=paste0("PC1 (",pca_data_perc[1],")"), 
       y=paste0("PC2 (",pca_data_perc[2],")"))



# color by condition/group
ggplot(df_pca_data, aes(PC1,PC2, color = condition))+
  geom_point(size=6)+
  labs(x=paste0("PC1 (",pca_data_perc[1],")"), 
       y=paste0("PC2 (",pca_data_perc[2],")"))+
  geom_text_repel(aes(label=sample),
                  point.padding = 20)


###############################################################################################################
##################################  DESEQ2 analysis ####################################################
###############################################################################################################


# R version 4.2.1 (2022-06-23)
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")  # v. 1.30.18
if (!require("DESeq2", quietly = TRUE))
  BiocManager::install("DESeq2")  # v. 1.36.0
library(DESeq2) # D.E.G.
FC <- 1.5 # Fold-change cutoff
FDR <- 0.1 # FDR cutoff
alpha <- 0.1 # independent filtering, default

#  Prepare data --------------------
raw_counts = read.csv("converted_counts_data.csv")
row.names(raw_counts) <- raw_counts$User_ID
raw_counts <- raw_counts[, -(1:3)] # delete 3 columns of IDs
str(raw_counts)

# Factors coded: Treatment --> A
col_data <- data.frame(
  "A" = c("447CONTROL", "447CONTROL", "447CONTROL", "447BACT", "447BACT", "447BACT", "447DOAB", "447DOAB", "447DOAB", "450CONTROL", "450CONTROL", "450CONTROL", "450BACT", "450BACT", "450BACT", "450DOAB", "450DOAB", "450DOAB")
)
row.names(col_data) <- colnames(raw_counts)
col_data

#Set reference level 
col_data[, 1] <- as.factor(col_data[, 1])
col_data[, 1] <- relevel(col_data[, 1], "447CONTROL")

# Run DESeq2--------------------
dds <- DESeq2::DESeqDataSetFromMatrix(
   countData = raw_counts,
   colData = col_data,
   design = ~  A 
)
dds = DESeq2::DESeq(dds) 

# Perform contrasts for paired comparisons of interest--------------------

# Comparison 1 of 6:  447BACT-447CONTROL
res <- DESeq2::results(dds,
  contrast = c("A", "447BACT", "447CONTROL"),
  lfcThreshold = log2(FC),
  altHypothesis = "greaterAbs",
  independentFiltering = TRUE,
  alpha = alpha
)
# Results
summary(res)
plotMA(res)
plotCounts(dds, gene = which.min(res$padj), intgroup = colnames(col_data)[2])
res <- subset(res, padj < FDR & abs(log2FoldChange) > log2(FC)) # Select
table(sign(res$log2FoldChange)) # N. of genes Down, Up
res <- res[order(-res$log2FoldChange), ] #sort
head(res) #top upregulated
tail(res) #top downregulated

# Comparison 2 of 6:  447DOAB-447CONTROL
res <- DESeq2::results(dds,
  contrast = c("A", "447DOAB", "447CONTROL"),
  lfcThreshold = log2(FC),
  altHypothesis = "greaterAbs",
  independentFiltering = TRUE,
  alpha = alpha
)
# Results
summary(res)
plotMA(res)
plotCounts(dds, gene = which.min(res$padj), intgroup = colnames(col_data)[2])
res <- subset(res, padj < FDR & abs(log2FoldChange) > log2(FC)) # Select
table(sign(res$log2FoldChange)) # N. of genes Down, Up
res <- res[order(-res$log2FoldChange), ] #sort
head(res) #top upregulated
tail(res) #top downregulated

# Comparison 3 of 6:  450BACT-447BACT
res <- DESeq2::results(dds,
  contrast = c("A", "450BACT", "447BACT"),
  lfcThreshold = log2(FC),
  altHypothesis = "greaterAbs",
  independentFiltering = TRUE,
  alpha = alpha
)
#Results 
summary(res)
plotMA(res)
plotCounts(dds, gene = which.min(res$padj), intgroup = colnames(col_data)[2])
res <- subset(res, padj < FDR & abs(log2FoldChange) > log2(FC)) # Select
table(sign(res$log2FoldChange)) # N. of genes Down, Up
res <- res[order(-res$log2FoldChange), ] #sort
head(res) #top upregulated
tail(res) #top downregulated

# Comparison 4 of 6:  450BACT-450CONTROL
res <- DESeq2::results(dds,
  contrast = c("A", "450BACT", "450CONTROL"),
  lfcThreshold = log2(FC),
  altHypothesis = "greaterAbs",
  independentFiltering = TRUE,
  alpha = alpha
)
# Results 
summary(res)
plotMA(res)
plotCounts(dds, gene = which.min(res$padj), intgroup = colnames(col_data)[2])
res <- subset(res, padj < FDR & abs(log2FoldChange) > log2(FC)) # Select
table(sign(res$log2FoldChange)) # N. of genes Down, Up
res <- res[order(-res$log2FoldChange), ] #sort
head(res) #top upregulated
tail(res) #top downregulated

# Comparison 5 of 6:  450DOAB-447DOAB
res <- DESeq2::results(dds,
  contrast = c("A", "450DOAB", "447DOAB"),
  lfcThreshold = log2(FC),
  altHypothesis = "greaterAbs",
  independentFiltering = TRUE,
  alpha = alpha
)
# Results
summary(res)
plotMA(res)
plotCounts(dds, gene = which.min(res$padj), intgroup = colnames(col_data)[2])
res <- subset(res, padj < FDR & abs(log2FoldChange) > log2(FC)) # Select
table(sign(res$log2FoldChange)) # N. of genes Down, Up
res <- res[order(-res$log2FoldChange), ] #sort
head(res) #top upregulated
tail(res) #top downregulated

# Comparison 6 of 6:  450DOAB-450CONTROL
res <- DESeq2::results(dds,
  contrast = c("A", "450DOAB", "450CONTROL"),
  lfcThreshold = log2(FC),
  altHypothesis = "greaterAbs",
  independentFiltering = TRUE,
  alpha = alpha
)
# Results 
summary(res)
plotMA(res)
plotCounts(dds, gene = which.min(res$padj), intgroup = colnames(col_data)[2])
res <- subset(res, padj < FDR & abs(log2FoldChange) > log2(FC)) # Select
table(sign(res$log2FoldChange)) # N. of genes Down, Up
res <- res[order(-res$log2FoldChange), ] #sort
head(res) #top upregulated
tail(res) #top downregulated


###############################################################################################################
##################################   ENHANCED VOLCANO PLOTS ####################################################
###############################################################################################################

if (!requireNamespace('BiocManager', quietly = TRUE))
  install.packages('BiocManager')
BiocManager::install('EnhancedVolcano')
install.packages("rlang")
install.packages("ggrepel")
library(ggplot2)
library(EnhancedVolcano)

r6 <- read.csv("447BACTvs447Ctl.csv", header=T)
r6 <- r6 %>% 
  column_to_rownames(var = "NCBIID")

r5 <- read.csv("447DOABvs447Ctl.csv", header=T)
r5 <- r5 %>% 
  column_to_rownames(var = "NCBIID")
              
r4 <- read.csv("450DOABvs447DOAB.csv", header=T)
r4 <- r4 %>% 
  column_to_rownames(var = "NCBIID")


r3 <- read.csv("450BACTvs447BACT.csv", header=T)
r3 <- r3 %>% 
  column_to_rownames(var = "NCBIID")


r2 <- read.csv("450DOABvs450CTL.csv", header=T)
r2 <- r2 %>% 
  column_to_rownames(var = "NCBIID")


r1 <- read.csv("450BACTvs450CTL.csv", header=T)
r1 <- r1 %>% 
  column_to_rownames(var = "NCBIID")



p1 <- EnhancedVolcano(r1,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'AdjPval',
                      title ="450BACT vs 450CTL",
                      titleLabSize = 12,
                      subtitle = "",
                      pCutoff = 10e-2,
                      FCcutoff=1.5,
                      xlim = c(-8, 8),
                      pointSize = 1,
                      axisLabSize = 14,
                      gridlines.major=TRUE,
                      gridlines.minor=FALSE,
                      border = "full",
                      borderWidth = 0.5,
                      caption = "",
                      legendPosition = '')

p2 <- EnhancedVolcano(r2,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'AdjPval',
                      title ="450DOAB vs 450CTL",
                      titleLabSize = 12,
                      subtitle = "",
                      pCutoff = 10e-2,
                      FCcutoff=1.5,
                      xlim = c(-8, 8),
                      pointSize = 1,
                      axisLabSize = 14,
                      gridlines.major=TRUE,
                      gridlines.minor=FALSE,
                      border = "full",
                      borderWidth = 0.5,
                      caption = "",
                      legendPosition = '')
p3 <- EnhancedVolcano(r3,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'AdjPval',
                      title ="450BACT vs 447BACT",
                      titleLabSize = 12,
                      subtitle = "",
                      pCutoff = 10e-2,
                      FCcutoff=1.5,
                      xlim = c(-8, 8),
                      pointSize = 1,
                      axisLabSize = 14,
                      gridlines.major=TRUE,
                      gridlines.minor=FALSE,
                      border = "full",
                      borderWidth = 0.5,
                      caption = "",
                      legendPosition = '')
p4 <- EnhancedVolcano(r4,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'AdjPval',
                      title ="450DOAB vs 447DOAB",
                      titleLabSize = 12,
                      subtitle = "",
                      pCutoff = 10e-2,
                      FCcutoff=1.5,
                      xlim = c(-8, 8),
                      pointSize = 1,
                      axisLabSize = 14,
                      gridlines.major=TRUE,
                      gridlines.minor=FALSE,
                      border = "full",
                      borderWidth = 0.5,
                      caption = "",
                      legendPosition = '')
p5 <- EnhancedVolcano(r5,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'AdjPval',
                      title ="447DOAB vs 447CTL",
                      titleLabSize = 12,
                      subtitle = "",
                      pCutoff = 10e-2,
                      FCcutoff=1.5,
                      xlim = c(-8, 8),
                      pointSize = 1,
                      axisLabSize = 14,
                      gridlines.major=TRUE,
                      gridlines.minor=FALSE,
                      border = "full",
                      borderWidth = 0.5,
                      caption = "",
                      legendPosition = '')
p6 <- EnhancedVolcano(r6,
                      lab = NA,
                      x = 'log2FoldChange',
                      y = 'AdjPval',
                      title = "447BACT vs 447CTL",
                      titleLabSize = 12,
                      subtitle = "",
                      pCutoff = 10e-2,
                      FCcutoff=1.5,
                      xlim = c(-8, 8),
                      pointSize = 1,
                      axisLabSize = 14,
                      gridlines.major=TRUE,
                      gridlines.minor=FALSE,
                      border = "full",
                      borderWidth = 0.5,
                      caption = "",
                      legendPosition = '')

library(gridExtra)

grid.arrange(p3, p4, p5, p6, ncol = 2, nrow = 2)
grid.arrange(p1, p2, ncol = 2, nrow = 1)

