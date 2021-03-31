# Loading of required packages ----
# (Install the packages if not already done)
library(limma)
library(edgeR)
library(ggplot2)
library(DESeq2)
library(VennDiagram)
library(gplots)
library(spatstat)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(fgsea)
library(clusterProfiler)
library(biomaRt)

# Set Working Directory ----
setwd("Y:/Bernard/Biomechanics research unit/Research/R Working directory/Lymphangio/Data_and_Codes")
list.files()

# Upload data ----

# Read file
data.df <- read.table("RawCounts_all_LEC_II4_VEGFC_24h.txt", header=TRUE, sep="\t", na.strings="NA",
                      row.names=1)
str(data.df)
head(data.df)
tail(data.df)
nrow(data.df)

# Are there missing data ?
missing <- is.na(data.df)
head(missing)
sum(missing[,]==TRUE)
# => No missing data

# Rename data frame and columns
data_LEC_CTRL_TUM_VEGFC <- data.df
colnames(data_LEC_CTRL_TUM_VEGFC) <- c("LEC_CTRL_1", "LEC_CTRL_2", "LEC_CTRL_3", 
                                       "LEC_TUM_1", "LEC_TUM_2", "LEC_TUM_3", "LEC_TUM_4",
                                       "LEC_VEGFC_1", "LEC_VEGFC_2", "LEC_VEGFC_3", "LEC_VEGFC_4")

head(data_LEC_CTRL_TUM_VEGFC)
nrow(data_LEC_CTRL_TUM_VEGFC)
# There are 26334 gene IDs (mapped nucleotidic sequences on the reference genome)
# Save the table
write.table(data_LEC_CTRL_TUM_VEGFC, "dataframe_rawcounts_LEC_CTRL_TUM_VEGFC.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)


# Data exploration ----

# Total number of reads per sample
attach(data_LEC_CTRL_TUM_VEGFC)
TotalNumbReads <- c(sum(LEC_CTRL_1), sum(LEC_CTRL_2), sum(LEC_CTRL_3), 
                    sum(LEC_TUM_1), sum(LEC_TUM_2), sum(LEC_TUM_3), sum(LEC_TUM_4), 
                    sum(LEC_VEGFC_1), sum(LEC_VEGFC_2), sum(LEC_VEGFC_3), sum(LEC_VEGFC_4))
TotalNumbReads
detach(data_LEC_CTRL_TUM_VEGFC)
# => The number of reads is not the same across samples. 
# => We will need to carry out normalization (scaling factor)

# Make a table of pseudocounts
pseudocounts <- log2(data_LEC_CTRL_TUM_VEGFC + 1)

# MA-plots between samples
x = pseudocounts[, 1]
y = pseudocounts[, 2]
# Change column numbers to check for the different comparisons between samples
# M-values
M = x - y
# A-values
A = (x + y)/2
df = data.frame(A, M)
ggplot(df, aes(x = A, y = M)) + geom_point(size = 1.5, alpha = 0.3) +
  geom_hline(yintercept = 0, color = "blue3") + 
  stat_smooth(se = FALSE, method = "loess", color = "red3")

# PCA (using DESeq2 package)
countData <- as.matrix(data_LEC_CTRL_TUM_VEGFC)
condition <- factor(c("control", "control", "control", "TUM", "TUM", "TUM", "TUM", 
                      "VEGFC", "VEGFC", "VEGFC", "VEGFC"))
dds <- DESeqDataSetFromMatrix(countData, DataFrame(condition), ~ condition)
rld <- rlog(dds)
plotPCA(rld)

# MDS plot
condition <- factor(c("control", "control", "control", "TUM", "TUM", "TUM", "TUM", 
                      "VEGFC", "VEGFC", "VEGFC", "VEGFC"))
plotMDS(pseudocounts, labels = condition)


# DIFFERENTIAL EXPRESSION ANALYSIS ----
# Using edgeR package (parametric approach)

# Experimental design ----
# There are three groups : control &  tumoral secretome (TUM) & VEGF-C (VEGFC).
# There are 4 biological replicates for TUM & VEGF and 3 biological replicates in the control group
group <- factor(c(1, 1, 1, 2, 2, 2, 2, 3, 3, 3, 3))
group
design <- model.matrix(~group)
design

# Filtering ----  
#  = filter out lowly expressed transcripts ahead of count-based differential expression analysis. 
# First create a DGEList
data_LEC_CTRL_TUM_VEGFC_DGEList <- DGEList(counts=data_LEC_CTRL_TUM_VEGFC, group=group)
data_LEC_CTRL_TUM_VEGFC_DGEList$samples

# Method used for filtering: 
# using counts per million (cpm) : treshold chosen at cpm = 0.25
# meaning a minimum of 4 reads in the smallest library 
# and a minimum of 5 reads in the biggest library
# in at least three of the samples
keep <- rowSums(cpm(data_LEC_CTRL_TUM_VEGFC_DGEList) > 0.25) >= 3
head(keep)
data_LEC_CTRL_TUM_VEGFC_Filtered <- data_LEC_CTRL_TUM_VEGFC_DGEList[keep, , keep.lib.sizes=FALSE]
head(data_LEC_CTRL_TUM_VEGFC_Filtered)
tail(data_LEC_CTRL_TUM_VEGFC_Filtered)
nrow(data_LEC_CTRL_TUM_VEGFC_Filtered)
# We keep 14330 rows (genes)

# Normalization ----
# in edgeR, normalization method is TMM.
data_LEC_CTRL_TUM_VEGFC_Norm <- calcNormFactors(data_LEC_CTRL_TUM_VEGFC_Filtered)
data_LEC_CTRL_TUM_VEGFC_Norm$samples

# Dispersion estimation ----
group
design
z <- estimateDisp(data_LEC_CTRL_TUM_VEGFC_Norm, design)
z$common.dispersion
plotBCV(z)

# Differential expression testing ----
fit <- glmQLFit(z, design)
# To compare 2 (TUM) versus 1 (CTRL)
qlf.2vs1 <- glmQLFTest(fit, coef=2)
topTags(qlf.2vs1, adjust.method="BH", sort.by="PValue")
# To compare 3 (VEGF) versus 1 (CTRL)
qlf.3vs1 <- glmQLFTest(fit, coef=3)
topTags(qlf.3vs1)
# To compare 3 (VEGF) versus 2 (TUM)
qlf.3vs2 <- glmQLFTest(fit, contrast=c(0,-1,1))
topTags(qlf.3vs2)

# MA plots between conditions ----
# MAplot CTRL(naive LEC) - TUM (teLEC)
de <- decideTestsDGE(qlf.2vs1, p.value = 0.05)
de
summary(de)
de.genes <- rownames(qlf.2vs1)[as.logical(de)]
plotSmear(qlf.2vs1, de.tags = de.genes, cex = 0.5, 
          main = "naive LEC vs teLEC",
          xlab = "Average expression (logCPM)",
          ylab= "log2 fold change",
          ylim= c(-10, 10))
abline(h = c(-1, 1), col = "blue")

# MAplot CTRL (naive LEC) - VEGFC (VEGF-LEC)
de <- decideTestsDGE(qlf.3vs1, p.value = 0.05)
de
summary(de)
de.genes <- rownames(qlf.3vs1)[as.logical(de)]
plotSmear(qlf.3vs1, de.tags = de.genes, cex = 0.5, 
          main = "naive LEC vs VEGF-LEC",
          xlab = "Average expression (logCPM)",
          ylab= "log2 fold change",
          ylim= c(-10, 10))
abline(h = c(-1, 1), col = "blue")

# MAplot TUM (teLEC) - VEGFC (VEGF-LEC)
de <- decideTestsDGE(qlf.3vs2, p.value = 0.05)
de
summary(de)
de.genes <- rownames(qlf.3vs2)[as.logical(de)]
plotSmear(qlf.3vs2, de.tags = de.genes, cex = 0.5, 
          main = "teLEC vs VEGF-LEC",
          xlab = "Average expression (logCPM)",
          ylab= "log2 fold change",
          ylim= c(-10, 10))
abline(h = c(-1, 1), col = "blue")

# Adding adjusted p-values to dataframes ----
# (correction for multiple comparisons using Benjamini & Hochberg method)
# Adjusted p-values for TUM versus CTRL
genes <- as.data.frame(qlf.2vs1$table)
head(genes)
FDR <- as.data.frame(p.adjust(qlf.2vs1$table$PValue, method="BH"))
str(FDR)
head(FDR)
colnames(FDR) <- "padj"
str(FDR)
str(genes)
DE_LEC_CTRL_TUM <- data.frame(genes, FDR)
head(DE_LEC_CTRL_TUM)
str(DE_LEC_CTRL_TUM)

# Adjusted p-values for VEGFC versus CTRL
genes <- as.data.frame(qlf.3vs1$table)
head(genes)
FDR <- as.data.frame(p.adjust(qlf.3vs1$table$PValue, method="BH"))
str(FDR)
head(FDR)
colnames(FDR) <- "padj"
str(FDR)
str(genes)
DE_LEC_CTRL_VEGFC <- data.frame(genes, FDR)
head(DE_LEC_CTRL_VEGFC)
str(DE_LEC_CTRL_VEGFC)

# Selecting up and down DE genes and save output tables ----
# naive LEC (CTRL) vs teLEC (TUM)
# downregulated and upregulated genes with logFC >= 1 or <= -1 and padj < 0.05
DE_LEC_CTRL_TUM_down <- DE_LEC_CTRL_TUM[DE_LEC_CTRL_TUM$logFC<=-1 & DE_LEC_CTRL_TUM$padj < 0.05,]
DE_LEC_CTRL_TUM_up <- DE_LEC_CTRL_TUM[DE_LEC_CTRL_TUM$logFC>=1 & DE_LEC_CTRL_TUM$padj < 0.05,]
DE_LEC_CTRL_TUM_up_down <- DE_LEC_CTRL_TUM[abs(DE_LEC_CTRL_TUM$logFC)>=1 & DE_LEC_CTRL_TUM$padj < 0.05,]
nrow(DE_LEC_CTRL_TUM_down)
# 539 genes
nrow(DE_LEC_CTRL_TUM_up)
# 520 genes
nrow(DE_LEC_CTRL_TUM_up_down)
# 1059 genes
# Save the tables with DE genes
write.table(DE_LEC_CTRL_TUM_down, "DE_LEC_CTRL_TUM_down.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(DE_LEC_CTRL_TUM_up, "DE_LEC_CTRL_TUM_up.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(DE_LEC_CTRL_TUM_up_down, "DE_LEC_CTRL_TUM_up_down.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(DE_LEC_CTRL_TUM, "DE_LEC_CTRL_TUM_all_genes.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

# naive LEC (CTRL) vs VEGF-LEC (VEGFC)
# downregulated and upregulated genes with logFC >= 1 or <= -1 and padj < 0.05
DE_LEC_CTRL_VEGFC_down <- DE_LEC_CTRL_VEGFC[DE_LEC_CTRL_VEGFC$logFC<=-1 & DE_LEC_CTRL_VEGFC$padj < 0.05,]
DE_LEC_CTRL_VEGFC_up <- DE_LEC_CTRL_VEGFC[DE_LEC_CTRL_VEGFC$logFC>=1 & DE_LEC_CTRL_VEGFC$padj < 0.05,]
DE_LEC_CTRL_VEGFC_up_down <- DE_LEC_CTRL_VEGFC[abs(DE_LEC_CTRL_VEGFC$logFC)>=1 & DE_LEC_CTRL_VEGFC$padj < 0.05,]
nrow(DE_LEC_CTRL_VEGFC_down)
# 777
nrow(DE_LEC_CTRL_VEGFC_up)
# 873
nrow(DE_LEC_CTRL_VEGFC_up_down)
# 1650
# Save the tables with DE genes
write.table(DE_LEC_CTRL_VEGFC_down, "DE_LEC_CTRL_VEGFC_down.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(DE_LEC_CTRL_VEGFC_up, "DE_LEC_CTRL_VEGFC_up.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(DE_LEC_CTRL_VEGFC_up_down, "DE_LEC_CTRL_VEGFC_up_down.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)
write.table(DE_LEC_CTRL_VEGFC, "DE_LEC_CTRL_VEGFC_all_genes.txt", sep="\t", row.names=TRUE, col.names=TRUE, quote=FALSE)

# Volcano plots ----

# naive LEC (CTRL) vs teLEC (TUM)
head(DE_LEC_CTRL_TUM)
# Define coordinates for up and down DE genes
gray <- subset(DE_LEC_CTRL_TUM, padj < 0.05)
blue <- subset(gray, logFC <= -1.0)
red <- subset(gray, logFC >= 1.0) 
gray <- subset(gray, abs(logFC)<1)
black <- subset(DE_LEC_CTRL_TUM, padj >= 0.05)
ygray <- -log10(gray$padj)
xgray <- gray$logFC
yblue <- -log10(blue$padj)
xblue <- blue$logFC
yred <- -log10(red$padj)
xred <- red$logFC
yblack <- -log10(black$padj)
xblack <- black$logFC
# Volcano plot
# blue = DE downregulated genes / red = DE upregulated genes 
# grey = DE genes with FDR controlled at 0.05 / black = equally expressed genes
plot(xgray, ygray, pch=16, cex=0.5, col="gray60",
     main="naive LEC vs teLEC",
     xlab= expression(paste(log[2],"fold change")), ylab= expression(paste(-log[10],"BH padj value")),
     xlim=c(-9, 9), ylim=c(0, 17))
points(xblack, yblack, pch=16, cex=0.5, col="black")
points(xblue, yblue, pch=16, cex=0.5, col="blue")
points(xred, yred, pch=16, cex=0.5, col="red")
abline(h=c(1.3), lwd=1.5, lty=2, col="gray60")
abline(v=c(-1, 1), lwd=1.5, lty=2, col= c("blue", "red"))
text(-8, 2, "FDR at 0.05",col="gray60", cex=0.8)
text(-6, 15, "539 down-regulated genes", col = "blue", cex =0.9)
text(6, 15, "520 up-regulated genes", col = "red", cex =0.9)

dev.off()

# naive LEC (CTRL) vs VEGF-LEC (VEGFC)
head(DE_LEC_CTRL_VEGFC)
# Define coordinates for up and down DE genes
gray <- subset(DE_LEC_CTRL_VEGFC, padj < 0.05)
blue <- subset(gray, logFC <= -1.0)
red <- subset(gray, logFC >= 1.0) 
gray <- subset(gray, abs(logFC)<1)
black <- subset(DE_LEC_CTRL_VEGFC, padj >= 0.05)
ygray <- -log10(gray$padj)
xgray <- gray$logFC
yblue <- -log10(blue$padj)
xblue <- blue$logFC
yred <- -log10(red$padj)
xred <- red$logFC
yblack <- -log10(black$padj)
xblack <- black$logFC
# Volcano plot
# blue = DE downregulated genes / red = DE upregulated genes 
# grey = DE genes with FDR controlled at 0.05 / black = equally expressed genes
plot(xgray, ygray, pch=16, cex=0.5, col="gray60",
     main="naive LEC vs VEGF-LEC",
     xlab= expression(paste(log[2],"fold change")), ylab= expression(paste(-log[10],"BH padj value")),
     xlim=c(-9, 9), ylim=c(0, 17))
points(xblack, yblack, pch=16, cex=0.5, col="black")
points(xblue, yblue, pch=16, cex=0.5, col="blue")
points(xred, yred, pch=16, cex=0.5, col="red")
abline(h=c(1.3), lwd=1.5, lty=2, col="gray60")
abline(v=c(-1, 1), lwd=1.5, lty=2, col= c("blue", "red"))
text(-8, 2, "FDR at 0.05",col="gray60", cex=0.8)
text(-6, 15, "777 down-regulated genes", col = "blue", cex =0.9)
text(6, 15, "873 up-regulated genes", col = "red", cex =0.9)

dev.off()

# POST-DIFFERENTIAL EXPRESSION ANALYSES ----
## Reload packages (if needed)
library(limma)
library(edgeR)
library(ggplot2)
library(VennDiagram)
library(gplots)
library(spatstat)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(fgsea)
library(clusterProfiler)
library(biomaRt)

## Set Working Directory (if needed)
setwd("Y:/Bernard/Biomechanics research unit/Research/R Working directory/Lymphangio/Data_and_Codes")
list.files()

# Venn diagrams comparing DE genes in teLEC and VEGF-LEC ----
# Use previously saved files
DE_LEC_CTRL_TUM_up <- read.table("DE_LEC_CTRL_TUM_up.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_TUM_up)
nrow(DE_LEC_CTRL_TUM_up)
DE_LEC_CTRL_TUM_down <- read.table("DE_LEC_CTRL_TUM_down.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_TUM_down)
nrow(DE_LEC_CTRL_TUM_down)
DE_LEC_CTRL_VEGFC_up <- read.table("DE_LEC_CTRL_VEGFC_up.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_VEGFC_up)
nrow(DE_LEC_CTRL_VEGFC_up)
DE_LEC_CTRL_VEGFC_down <- read.table("DE_LEC_CTRL_VEGFC_down.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_VEGFC_down)
nrow(DE_LEC_CTRL_VEGFC_down)

# Convert gene IDs as factors
DE_LEC_CTRL_TUM_up_id <- as.factor(rownames(DE_LEC_CTRL_TUM_up))
head(DE_LEC_CTRL_TUM_up_id)
DE_LEC_CTRL_TUM_down_id <- as.factor(rownames(DE_LEC_CTRL_TUM_down))
head(DE_LEC_CTRL_TUM_down_id)
DE_LEC_CTRL_VEGFC_up_id <- as.factor(rownames(DE_LEC_CTRL_VEGFC_up))
head(DE_LEC_CTRL_VEGFC_up_id)
DE_LEC_CTRL_VEGFC_down_id <- as.factor(rownames(DE_LEC_CTRL_VEGFC_down))
head(DE_LEC_CTRL_VEGFC_down_id)

# Calculate the number of common genes between conditions 
# between upregulated genes
length(intersect(DE_LEC_CTRL_TUM_up_id, DE_LEC_CTRL_VEGFC_up_id))
# between downregulated genes
length(intersect(DE_LEC_CTRL_TUM_down_id, DE_LEC_CTRL_VEGFC_down_id))

# Venn diagram with upregulated genes
venn.plot <- draw.pairwise.venn(520, 873, 159, category = c("teLEC", "VEGF-LEC"),
                                euler.d = FALSE, scaled = FALSE, inverted = FALSE, rotation.degree = 180,
                                col = c("gold", "palegreen3"), 
                                fill = c("gold", "palegreen3"), 
                                cex = rep(2, 3), cat.pos = c(0,0),
                                cat.cex = c(1.5, 1.5),
                                label.col = c("gold4", "darkolivegreen", "palegreen4"),
                                cat.col = c("gold4", "palegreen4"))
dev.off()

# Venn diagram with downregulated genes
venn.plot <- draw.pairwise.venn(539, 777, 243, category = c("teLEC", "VEGF-LEC"),
                                euler.d = FALSE, scaled = FALSE, inverted = FALSE, rotation.degree = 180,
                                col = c("gold", "palegreen3"), 
                                fill = c("gold", "palegreen3"), 
                                cex = rep(2, 3), cat.pos = c(0,0),
                                cat.cex = c(1.5, 1.5),
                                label.col = c("gold4", "darkolivegreen", "palegreen4"),
                                cat.col = c("gold4", "palegreen4"))
dev.off()

# Heatmap top 300 genes in each dataset (teLEC + VEGF-LEC) ----
# Top 300 of the genes ranked according to the absolute values of log2FC
# Representing scaled raw counts for each sample

# Take DE genes (up and down)
DE_LEC_CTRL_TUM_up_down <- read.table("DE_LEC_CTRL_TUM_up_down.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_TUM_up_down)
nrow(DE_LEC_CTRL_TUM_up_down)
DE_LEC_CTRL_VEGFC_up_down <- read.table("DE_LEC_CTRL_VEGFC_up_down.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_VEGFC_up_down)
nrow(DE_LEC_CTRL_VEGFC_up_down)

# Order by absolute fold change (and padj)
DE_LEC_CTRL_TUM_up_down_by_abs_logFC <- DE_LEC_CTRL_TUM_up_down[order(-abs(DE_LEC_CTRL_TUM_up_down$logFC), DE_LEC_CTRL_TUM_up_down$padj), ]
DE_LEC_CTRL_VEGFC_up_down_by_abs_logFC <- DE_LEC_CTRL_VEGFC_up_down[order(-abs(DE_LEC_CTRL_VEGFC_up_down$logFC),DE_LEC_CTRL_VEGFC_up_down$padj), ]

# Select the first 300 rows
Top300_DE_LEC_CTRL_TUM <- head(DE_LEC_CTRL_TUM_up_down_by_abs_logFC, 300)
Top300_DE_LEC_CTRL_VEGFC <- head(DE_LEC_CTRL_VEGFC_up_down_by_abs_logFC, 300)

# Take gene IDs and convert as factors
Top300_DE_LEC_CTRL_TUM_id <- as.factor(rownames(Top300_DE_LEC_CTRL_TUM))
head(Top300_DE_LEC_CTRL_TUM_id)
Top300_DE_LEC_CTRL_VEGFC_id <- as.factor(rownames(Top300_DE_LEC_CTRL_VEGFC))
head(Top300_DE_LEC_CTRL_VEGFC_id)

# Group all the genes (and make a dataframe)
top_DE_genes <- union(Top300_DE_LEC_CTRL_TUM_id, Top300_DE_LEC_CTRL_VEGFC_id)
length(top_DE_genes)
length(intersect(Top300_DE_LEC_CTRL_TUM_id, Top300_DE_LEC_CTRL_VEGFC_id))
# 540 genes in total (as we have 60 genes in common)
df_top_DE_genes <- as.data.frame(top_DE_genes)
row.names(df_top_DE_genes) <- df_top_DE_genes$top_DE_genes

# Merge this dataframe (top DE genes) with the initial table of raw counts
data_LEC_CTRL_TUM_VEGFC <- read.table("dataframe_rawcounts_LEC_CTRL_TUM_VEGFC.txt", 
                                      header=TRUE, sep="\t", na.strings="NA",row.names=1)

merge1 <- merge(df_top_DE_genes, data_LEC_CTRL_TUM_VEGFC, by = "row.names")
row.names(merge1) <- merge1$Row.names
head(merge1)
merge2 <- merge1[ , c(3:13)]
head(merge2)
colnames(merge2) <- c("naive LEC 1", "naive LEC 2", "naive LEC 3", 
                      "teLEC 1", "teLEC 2", "teLEC 3", "teLEC 4", 
                      "VEGF-LEC 1", "VEGF-LEC 2", "VEGF-LEC 3", "VEGF-LEC 4" )
head(merge2)

# Make the heatmap
# Convert dataframe to a matrix
data_mat <- data.matrix(merge2)
# heatmap.2 package is doing clustering first, before scaling #
# Yet we want scaling before clustering
# So let's scale the data before visualizing
# We can only scale the columns --> so we also need to transpose the matrix
data_scaled <- t(scale(t(data_mat)))
# Create his own palette of colors
my_palette <- colorRampPalette(c("blue", "white", "red"))
# Heatmap
heatmap.2(data_scaled, col = my_palette(n=25),
          trace = "none",
          Rowv = TRUE,
          Colv = FALSE,
          dendrogram = "row",
          density.info = "none",
          key = TRUE,
          keysize = 1,
          margins = c(8, 3),
          cexCol = 1.2,
          labRow = FALSE,
          ylab = "DE genes",
          xlab = "individual samples")
dev.off()

# Heatmaps top60 DE genes (30UP and 30DOWN) in each dataset (teLEC + VEGF-LEC) ----
# 30 most up- and 30 most down-regulated genes (according to log2FC values)
# Representing log2FC of these "top60" genes for each condition (teLEC or VEGF-LEC)

# Take DE genes (up and down)
DE_LEC_CTRL_TUM_up_down <- read.table("DE_LEC_CTRL_TUM_up_down.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_TUM_up_down)
nrow(DE_LEC_CTRL_TUM_up_down)
DE_LEC_CTRL_VEGFC_up_down <- read.table("DE_LEC_CTRL_VEGFC_up_down.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_VEGFC_up_down)
nrow(DE_LEC_CTRL_VEGFC_up_down)

# Order by decreasing fold change (and padj)
DE_LEC_CTRL_TUM_up_down_by_logFC <- DE_LEC_CTRL_TUM_up_down[order(-DE_LEC_CTRL_TUM_up_down$logFC,DE_LEC_CTRL_TUM_up_down$padj), ]
DE_LEC_CTRL_VEGFC_up_down_by_logFC <- DE_LEC_CTRL_VEGFC_up_down[order(-DE_LEC_CTRL_VEGFC_up_down$logFC,DE_LEC_CTRL_VEGFC_up_down$padj), ]
head(DE_LEC_CTRL_TUM_up_down_by_logFC)
head(DE_LEC_CTRL_VEGFC_up_down_by_logFC)

# Select the 30 first and 30 last rows (30 most up- and down-regulated genes)
Top30_UP_LEC_CTRL_TUM <- head(DE_LEC_CTRL_TUM_up_down_by_logFC, 30)
Top30_DOWN_LEC_CTRL_TUM <- tail(DE_LEC_CTRL_TUM_up_down_by_logFC, 30)
Top30_UP_LEC_CTRL_VEGFC <- head(DE_LEC_CTRL_VEGFC_up_down_by_logFC, 30)
Top30_DOWN_LEC_CTRL_VEGFC <- tail(DE_LEC_CTRL_VEGFC_up_down_by_logFC, 30)

# Take gene IDs and convert as factors
Top30_UP_LEC_CTRL_TUM_id <- as.factor(rownames(Top30_UP_LEC_CTRL_TUM))
Top30_DOWN_LEC_CTRL_TUM_id <- as.factor(rownames(Top30_DOWN_LEC_CTRL_TUM))
Top30_UP_LEC_CTRL_VEGFC_id <- as.factor(rownames(Top30_UP_LEC_CTRL_VEGFC))
Top30_DOWN_LEC_CTRL_VEGFC_id <- as.factor(rownames(Top30_DOWN_LEC_CTRL_VEGFC))

# Regroup top 30 up- and down-regulated genes (= "top 60") for teLEC (and make dataframe)
Top60_UP_DOWN_LEC_CTRL_TUM_id <- union(Top30_UP_LEC_CTRL_TUM_id, Top30_DOWN_LEC_CTRL_TUM_id)
df_Top60_UP_DOWN_LEC_CTRL_TUM <- as.data.frame(Top60_UP_DOWN_LEC_CTRL_TUM_id)
row.names(df_Top60_UP_DOWN_LEC_CTRL_TUM) <- df_Top60_UP_DOWN_LEC_CTRL_TUM$Top60_UP_DOWN_LEC_CTRL_TUM_id
head(df_Top60_UP_DOWN_LEC_CTRL_TUM)

# Regroup top 30 up- and down-regulated genes (= "top 60") for VEGF-LEC
Top60_UP_DOWN_LEC_CTRL_VEGFC_id <- union(Top30_UP_LEC_CTRL_VEGFC_id, Top30_DOWN_LEC_CTRL_VEGFC_id)
df_Top60_UP_DOWN_LEC_CTRL_VEGFC <- as.data.frame(Top60_UP_DOWN_LEC_CTRL_VEGFC_id)
row.names(df_Top60_UP_DOWN_LEC_CTRL_VEGFC) <- df_Top60_UP_DOWN_LEC_CTRL_VEGFC$Top60_UP_DOWN_LEC_CTRL_VEGFC_id
head(df_Top60_UP_DOWN_LEC_CTRL_VEGFC)

# Heatmap showing Top 60 DE genes in teLEC
# Merge the dataframe (top 60 DE genes) with the tables obtained after DE (all genes, logFC, pval, padj)
DE_LEC_CTRL_TUM_all_genes <- read.table("DE_LEC_CTRL_TUM_all_genes.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_TUM_all_genes)
nrow(DE_LEC_CTRL_TUM_all_genes)
DE_LEC_CTRL_VEGFC_all_genes <- read.table("DE_LEC_CTRL_VEGFC_all_genes.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_VEGFC_all_genes)
nrow(DE_LEC_CTRL_VEGFC_all_genes)
merge1 <- merge(df_Top60_UP_DOWN_LEC_CTRL_TUM, DE_LEC_CTRL_TUM_all_genes, by = "row.names", all.x = TRUE)
row.names(merge1) <- merge1$Row.names
head(merge1)
# Keep only logFC and padj values
merge1 <- merge1[ , c(3,7)]
head(merge1)
colnames(merge1) <- c("logFC_teLEC", "padj_teLEC")
head(merge1)
# Add results obtained for those genes in the VEGF-LEC condition (as comparison)
merge2 <- merge(merge1, DE_LEC_CTRL_VEGFC_all_genes, by = "row.names", all.x = TRUE)
row.names(merge2) <- merge2$Row.names
head(merge2)
# Keep only logFC and padj values
merge2 <- merge2[ , c(2:4,8)]
head(merge2)
colnames(merge2) <- c("logFC_teLEC", "padj_teLEC","logFC_VEGF-LEC", "padj_VEGF-LEC" )
head(merge2)
# If padj values > 0.05 (non-significance), set the logFC values to 0 (as these are non-DE genes)
x <- merge2
x$logFC_teLEC <- x$logFC_teLEC * (as.integer(x$padj_teLEC < 0.05))
x$logFC_teLEC
x$`logFC_VEGF-LEC` <- x$`logFC_VEGF-LEC` * (as.integer(x$`padj_VEGF-LEC` < 0.05))
x$`logFC_VEGF-LEC`
head(x)
# Keep only the log2FC
data_logFC <- x[ , c(1,3)]
head(data_logFC)
colnames(data_logFC) <- c("teLEC", "VEGF-LEC")
head(data_logFC)
# Order by decreasing log2FC in teLEC
data_logFC_ordered <- data_logFC[order(-data_logFC$teLEC), ]
head(data_logFC_ordered)
# Convert to a matrix
data_mat <- data.matrix(data_logFC_ordered)
# Make a vector with n+1 breaks
ncol <- 100
rampbreaks <- seq(-10, 10, length.out = ncol+1)
# Heatmap
pheatmap(data_mat, 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = rampbreaks,
         scale = "none",
         legend = TRUE,
         angle_col = 45,
         fontsize_row = 8,
         fontsize_col = 8,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         cellwidth = 40,
         treeheight_row = 0,
         display_numbers = TRUE)
dev.off()

# Heatmap showing Top 60 DE genes in VEGF-LEC
# Merge the dataframe (top 60 DE genes) with the tables obtained after DE (all genes, logFC, pval, padj)
DE_LEC_CTRL_TUM_all_genes <- read.table("DE_LEC_CTRL_TUM_all_genes.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_TUM_all_genes)
nrow(DE_LEC_CTRL_TUM_all_genes)
DE_LEC_CTRL_VEGFC_all_genes <- read.table("DE_LEC_CTRL_VEGFC_all_genes.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_VEGFC_all_genes)
nrow(DE_LEC_CTRL_VEGFC_all_genes)
# First merge/add results obtained for those genes in the teLEC condition (as comparison)
merge1 <- merge(df_Top60_UP_DOWN_LEC_CTRL_VEGFC, DE_LEC_CTRL_TUM_all_genes, by = "row.names", all.x = TRUE)
row.names(merge1) <- merge1$Row.names
head(merge1)
# Keep only logFC and padj values
merge1 <- merge1[ , c(3,7)]
head(merge1)
colnames(merge1) <- c("logFC_teLEC", "padj_teLEC")
head(merge1)
# Merge with the results in VEGF-LEC condition
merge2 <- merge(merge1, DE_LEC_CTRL_VEGFC_all_genes, by = "row.names", all.x = TRUE)
row.names(merge2) <- merge2$Row.names
head(merge2)
# Keep only logFC and padj values
merge2 <- merge2[ , c(2:4,8)]
head(merge2)
colnames(merge2) <- c("logFC_teLEC", "padj_teLEC","logFC_VEGF-LEC", "padj_VEGF-LEC" )
head(merge2)
# If padj values > 0.05 (non-significance), set the logFC values to 0 (as these are non-DE genes)
x <- merge2
x$logFC_teLEC <- x$logFC_teLEC * (as.integer(x$padj_teLEC < 0.05))
x$logFC_teLEC
x$`logFC_VEGF-LEC` <- x$`logFC_VEGF-LEC` * (as.integer(x$`padj_VEGF-LEC` < 0.05))
x$`logFC_VEGF-LEC`
head(x)
# Keep only the log2FC
data_logFC <- x[ , c(1,3)]
head(data_logFC)
colnames(data_logFC) <- c("teLEC", "VEGF-LEC")
head(data_logFC)
# Order by decreasing log2FC in VEGF-LEC
data_logFC_ordered <- data_logFC[order(-data_logFC$`VEGF-LEC`), ]
head(data_logFC_ordered)
# Convert to a matrix
data_mat <- data.matrix(data_logFC_ordered)
# Make a vector with n+1 breaks
ncol <- 100
rampbreaks <- seq(-10, 10, length.out = ncol+1)
# Heatmap
pheatmap(data_mat, 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = rampbreaks,
         scale = "none",
         legend = TRUE,
         angle_col = 45,
         fontsize_row = 8,
         fontsize_col = 8,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         cellwidth = 40,
         treeheight_row = 0,
         display_numbers = TRUE)
dev.off()

# DE genes modulated in an opposite way in teLEC and VEGF-LEC ----

# Take de genes (up and down)
DE_LEC_CTRL_TUM_up <- read.table("DE_LEC_CTRL_TUM_up.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_TUM_up)
nrow(DE_LEC_CTRL_TUM_up)
DE_LEC_CTRL_VEGFC_up <- read.table("DE_LEC_CTRL_VEGFC_up.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_VEGFC_up)
nrow(DE_LEC_CTRL_VEGFC_up)
DE_LEC_CTRL_TUM_down <- read.table("DE_LEC_CTRL_TUM_down.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_TUM_down)
nrow(DE_LEC_CTRL_TUM_down)
DE_LEC_CTRL_VEGFC_down <- read.table("DE_LEC_CTRL_VEGFC_down.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_VEGFC_down)
nrow(DE_LEC_CTRL_VEGFC_down)

# Take gene IDs
DE_LEC_CTRL_TUM_up_id <- as.factor(rownames(DE_LEC_CTRL_TUM_up))
head(DE_LEC_CTRL_TUM_up_id)
DE_LEC_CTRL_VEGFC_up_id <- as.factor(rownames(DE_LEC_CTRL_VEGFC_up))
head(DE_LEC_CTRL_VEGFC_up_id)
DE_LEC_CTRL_TUM_down_id <- as.factor(rownames(DE_LEC_CTRL_TUM_down))
head(DE_LEC_CTRL_TUM_down_id)
DE_LEC_CTRL_VEGFC_down_id <- as.factor(rownames(DE_LEC_CTRL_VEGFC_down))
head(DE_LEC_CTRL_VEGFC_down_id)

# Select genes downregulated in teLEC and upregulated in VEGF-LEC 
DE_LEC_TUM_down_VEGFC_up <- intersect(DE_LEC_CTRL_TUM_down_id, DE_LEC_CTRL_VEGFC_up_id)
length(DE_LEC_TUM_down_VEGFC_up) # 23 genes

# Select genes upregulated in teLEC and downregulated in VEGF-LEC
DE_LEC_TUM_up_VEGFC_down <- intersect(DE_LEC_CTRL_TUM_up_id, DE_LEC_CTRL_VEGFC_down_id)
length(DE_LEC_TUM_up_VEGFC_down) # 27 genes

# Regroup both sets of genes and create a dataframe
DE_LEC_CTRL_TUM_VEGFC_opposites <- union(DE_LEC_TUM_down_VEGFC_up, DE_LEC_TUM_up_VEGFC_down)
length(DE_LEC_CTRL_TUM_VEGFC_opposites) # 50 genes
DE_LEC_CTRL_TUM_VEGFC_opposites_df <- as.data.frame(DE_LEC_CTRL_TUM_VEGFC_opposites)
row.names(DE_LEC_CTRL_TUM_VEGFC_opposites_df) <- DE_LEC_CTRL_TUM_VEGFC_opposites_df$DE_LEC_CTRL_TUM_VEGFC_opposites
head(DE_LEC_CTRL_TUM_VEGFC_opposites_df)

# Merge the "opposites/contrasting genes" with DE tables  (logFC, pval, padj) for both teLEC and VEGF-LEC
DE_LEC_CTRL_TUM_up_down <- read.table("DE_LEC_CTRL_TUM_up_down.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_TUM_up_down)
DE_LEC_CTRL_VEGFC_up_down <- read.table("DE_LEC_CTRL_VEGFC_up_down.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_VEGFC_up_down)
merge1 <- merge(DE_LEC_CTRL_TUM_VEGFC_opposites_df, DE_LEC_CTRL_TUM_up_down, by="row.names", all.x = TRUE)
row.names(merge1) <- merge1$Row.names
merge1 <- merge1[ , c(3,7)]
head(merge1)
colnames(merge1) <- c("logFC_teLEC", "padj_teLEC")
merge2 <- merge(merge1, DE_LEC_CTRL_VEGFC_up_down, by = "row.names", all.x = TRUE)
row.names(merge2) <- merge2$Row.names
head(merge2)
merge2 <- merge2[ , c(2:4,8)]
colnames(merge2) <- c("logFC_teLEC", "padj_teLEC", "logFC_VEGF-LEC", "padj_VEGF-LEC")
# Order according to log2FC in teLEC
contrasting_genes <- merge2[order(-merge2$logFC_teLEC), ]
head(contrasting_genes)

# Prepare dotchart for teLEC part
# create "upteLEC" column --> TRUE if logFC_teLEC >=0 (FALSE if < 0)
teLEC_contrasting_genes <- contrasting_genes
teLEC_contrasting_genes$upteLEC <- teLEC_contrasting_genes$logFC_teLEC >= 0 
teLEC_contrasting_genes
# create group column -> if upteLEC = TRUE then "teLEC up", if upteLEC = FALSE then "teLEC down"
teLEC_contrasting_genes$group[teLEC_contrasting_genes$upteLEC == TRUE] <- "teLEC up"
teLEC_contrasting_genes$group[teLEC_contrasting_genes$upteLEC == FALSE] <- "teLEC down"
teLEC_contrasting_genes
# make group columns as a factor
teLEC_contrasting_genes$group <- factor(teLEC_contrasting_genes$group)
str(teLEC_contrasting_genes)
# create color column --> red for teLEC up and blue for teLEC down
teLEC_contrasting_genes$color[teLEC_contrasting_genes$group == "teLEC up"] <- "red"
teLEC_contrasting_genes$color[teLEC_contrasting_genes$group == "teLEC down"] <- "blue"
str(teLEC_contrasting_genes)
head(teLEC_contrasting_genes)

# Prepare dotchart for VEGF-LEC part
# create "upVEGF-LEC" column --> TRUE if logFC_VEGF-LEC >=0 (FALSE if < 0)
VEGFLEC_contrasting_genes <- contrasting_genes
VEGFLEC_contrasting_genes$upVEGFLEC <- VEGFLEC_contrasting_genes$`logFC_VEGF-LEC` >= 0 
VEGFLEC_contrasting_genes
# create group column -> if upVEGFLEC = TRUE then "VEGF-LEC up", if upVEGFLEC = FALSE then "VEGF-LEC down"
VEGFLEC_contrasting_genes$group[VEGFLEC_contrasting_genes$upVEGFLEC == TRUE] <- "VEGF-LEC up"
VEGFLEC_contrasting_genes$group[VEGFLEC_contrasting_genes$upVEGFLEC == FALSE] <- "VEGF-LEC down"
VEGFLEC_contrasting_genes
# make group columns as a factor
VEGFLEC_contrasting_genes$group <- factor(VEGFLEC_contrasting_genes$group)
str(VEGFLEC_contrasting_genes)
# reorder levels of factor (group column)
VEGFLEC_contrasting_genes$group <- relevel(VEGFLEC_contrasting_genes$group, "VEGF-LEC up")
str(VEGFLEC_contrasting_genes)
# create color column --> red for VEGF-LEC up and blue for VEGF-LEC down
VEGFLEC_contrasting_genes$color[VEGFLEC_contrasting_genes$group == "VEGF-LEC up"] <- "red"
VEGFLEC_contrasting_genes$color[VEGFLEC_contrasting_genes$group == "VEGF-LEC down"] <- "blue"
str(VEGFLEC_contrasting_genes)
head(VEGFLEC_contrasting_genes)

# make the two dotcharts
dotchart(	
  teLEC_contrasting_genes$logFC_teLEC,
  labels = row.names(teLEC_contrasting_genes), cex=0.6,
  groups = teLEC_contrasting_genes$group,
  gcolor = "black",
  color = teLEC_contrasting_genes$color,
  pch = 16,
  main = "teLEC", cex.main = 1.5,
  xlab = "log2FC in teLEC",
  cex.lab = 1.2,
  xlim = c(-5,5)
)
dev.off()

dotchart(	
  VEGFLEC_contrasting_genes$`logFC_VEGF-LEC`,
  labels = row.names(VEGFLEC_contrasting_genes), cex=0.6,
  groups = VEGFLEC_contrasting_genes$group,
  gcolor = "black",
  color = VEGFLEC_contrasting_genes$color,
  pch = 16,
  main = "VEGF-LEC", cex.main = 1.5,
  xlab = "log2FC in VEGF-LEC",
  cex.lab = 1.2,
  xlim = c(-5,5),
)
dev.off()

# Gene Set Enrichment Analysis ----
# Using fgsea package

# Use collection database (MSigDB, Broad Institute):
# download the .gmt file in the appropriate directory for R usage
# upload here the different files
# 1) GO Biological Process : c5.bp.v7.0.symbols.gmt
GOBP_file <- "Y:\\Bernard\\Biomechanics research unit\\Research\\R Working directory\\Niki\\c5.bp.v7.0.symbols.gmt"
GOBP <- gmtPathways(GOBP_file)
# 2) CP KEGG gene sets : c2.cp.kegg.v7.0.symbols.gmt
KEGG_file <- "Y:\\Bernard\\Biomechanics research unit\\Research\\R Working directory\\Niki\\c2.cp.kegg.v7.0.symbols.gmt" 
KEGG <- gmtPathways(KEGG_file)
# 3) CP Reactome : c2.cp.reactome.v7.1.symbols.gmt
REACTOME_file <- "Y:\\Bernard\\Biomechanics research unit\\Research\\R Working directory\\Niki\\c2.cp.reactome.v7.1.symbols.gmt"
REACTOME <- gmtPathways(REACTOME_file)

# Create gene lists for enrichment analysis with fgsea package
DE_LEC_CTRL_TUM_all_genes <- read.table("DE_LEC_CTRL_TUM_all_genes.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_TUM_all_genes)
DE_LEC_CTRL_VEGFC_all_genes <- read.table("DE_LEC_CTRL_VEGFC_all_genes.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_VEGFC_all_genes)
gene_list_LEC_CTRL_TUM = DE_LEC_CTRL_TUM_all_genes$logFC
gene_list_LEC_CTRL_VEGFC = DE_LEC_CTRL_VEGFC_all_genes$logFC
names(gene_list_LEC_CTRL_TUM) = row.names(DE_LEC_CTRL_TUM_all_genes)  
names(gene_list_LEC_CTRL_VEGFC) = row.names(DE_LEC_CTRL_VEGFC_all_genes) 
gene_list_LEC_CTRL_TUM = sort(gene_list_LEC_CTRL_TUM, decreasing = TRUE)
gene_list_LEC_CTRL_VEGFC = sort(gene_list_LEC_CTRL_VEGFC, decreasing = TRUE)
gene_list_LEC_CTRL_TUM = gene_list_LEC_CTRL_TUM[!duplicated(names(gene_list_LEC_CTRL_TUM))]
gene_list_LEC_CTRL_VEGFC = gene_list_LEC_CTRL_VEGFC[!duplicated(names(gene_list_LEC_CTRL_VEGFC))]
head(gene_list_LEC_CTRL_TUM, 50)
head(gene_list_LEC_CTRL_VEGFC, 50)
str(gene_list_LEC_CTRL_TUM)
str(gene_list_LEC_CTRL_VEGFC)

# GSEA with fgsea package
# with GO BP database
fgseaRes_GOBP_TUM <- fgseaMultilevel(pathways = GOBP, 
                                     stats = gene_list_LEC_CTRL_TUM, 
                                     minSize=15, 
                                     maxSize = 500)

fgseaRes_GOBP_VEGFC <- fgseaMultilevel(pathways = GOBP, 
                                       stats = gene_list_LEC_CTRL_VEGFC, 
                                       minSize=15, 
                                       maxSize = 500)

# with KEGG database
fgseaRes_KEGG_TUM <- fgseaMultilevel(pathways = KEGG, 
                                     stats = gene_list_LEC_CTRL_TUM, 
                                     minSize=15, 
                                     maxSize = 500)

fgseaRes_KEGG_VEGFC <- fgseaMultilevel(pathways = KEGG, 
                                       stats = gene_list_LEC_CTRL_VEGFC, 
                                       minSize=15, 
                                       maxSize = 500)
# with REACTOME database
fgseaRes_REACTOME_TUM <- fgseaMultilevel(pathways = REACTOME, 
                                         stats = gene_list_LEC_CTRL_TUM, 
                                         minSize=15, 
                                         maxSize = 500)

fgseaRes_REACTOME_VEGFC <- fgseaMultilevel(pathways = REACTOME, 
                                           stats = gene_list_LEC_CTRL_VEGFC, 
                                           minSize=15, 
                                           maxSize = 500)

# Select the 10 most upregulated pathways (positive enrichment score)
# Ranking based on padj and NES
topPathwaysUp_GOBP_TUM <- fgseaRes_GOBP_TUM[ES > 0][head(order(padj, -NES), n=10), pathway]
topPathwaysUp_GOBP_TUM
topPathwaysUp_GOBP_VEGFC <- fgseaRes_GOBP_VEGFC[ES > 0][head(order(padj, -NES), n=10), pathway]
topPathwaysUp_GOBP_VEGFC
topPathwaysUp_KEGG_TUM <- fgseaRes_KEGG_TUM[ES > 0][head(order(padj, -NES), n=10), pathway]
topPathwaysUp_KEGG_TUM
topPathwaysUp_KEGG_VEGFC <- fgseaRes_KEGG_VEGFC[ES > 0][head(order(padj, -NES), n=10), pathway]
topPathwaysUp_KEGG_VEGFC
topPathwaysUp_REACTOME_TUM <- fgseaRes_REACTOME_TUM[ES > 0][head(order(padj, -NES), n=10), pathway]
topPathwaysUp_REACTOME_TUM
topPathwaysUp_REACTOME_VEGFC <- fgseaRes_REACTOME_VEGFC[ES > 0][head(order(padj, -NES), n=10), pathway]
topPathwaysUp_REACTOME_VEGFC

# Select the 10 most downregulated pathways (negative enrichment score)
# Ranking based on padj and NES
topPathwaysDown_GOBP_TUM <- fgseaRes_GOBP_TUM[ES < 0][head(order(padj, NES), n=10), pathway]
topPathwaysDown_GOBP_TUM
topPathwaysDown_GOBP_VEGFC <- fgseaRes_GOBP_VEGFC[ES < 0][head(order(padj, NES), n=10), pathway]
topPathwaysDown_GOBP_VEGFC
topPathwaysDown_KEGG_TUM <- fgseaRes_KEGG_TUM[ES < 0][head(order(padj, NES), n=10), pathway]
topPathwaysDown_KEGG_TUM
topPathwaysDown_KEGG_VEGFC <- fgseaRes_KEGG_VEGFC[ES < 0][head(order(padj, NES), n=10), pathway]
topPathwaysDown_KEGG_VEGFC
topPathwaysDown_REACTOME_TUM <- fgseaRes_REACTOME_TUM[ES < 0][head(order(padj, NES), n=10), pathway]
topPathwaysDown_REACTOME_TUM
topPathwaysDown_REACTOME_VEGFC <- fgseaRes_REACTOME_VEGFC[ES < 0][head(order(padj, NES), n=10), pathway]
topPathwaysDown_REACTOME_VEGFC

# Plot the tables (10 most up- and down-regulated pathways using GOBP, KEGG or REACTOME)
plotGseaTable(GOBP[topPathwaysUp_GOBP_TUM], gene_list_LEC_CTRL_TUM, fgseaRes_GOBP_TUM, 
              gseaParam = 0.5)
dev.off()
plotGseaTable(GOBP[topPathwaysUp_GOBP_VEGFC], gene_list_LEC_CTRL_VEGFC, fgseaRes_GOBP_VEGFC, 
              gseaParam = 0.5)
dev.off()
plotGseaTable(KEGG[topPathwaysUp_KEGG_TUM], gene_list_LEC_CTRL_TUM, fgseaRes_KEGG_TUM, 
              gseaParam = 0.5)
dev.off()
plotGseaTable(KEGG[topPathwaysUp_KEGG_VEGFC], gene_list_LEC_CTRL_VEGFC, fgseaRes_KEGG_VEGFC, 
              gseaParam = 0.5)
dev.off()
plotGseaTable(REACTOME[topPathwaysUp_REACTOME_TUM], gene_list_LEC_CTRL_TUM, fgseaRes_REACTOME_TUM, 
              gseaParam = 0.5)
dev.off()
plotGseaTable(REACTOME[topPathwaysUp_REACTOME_VEGFC], gene_list_LEC_CTRL_VEGFC, fgseaRes_REACTOME_VEGFC, 
              gseaParam = 0.5)
dev.off()
plotGseaTable(GOBP[topPathwaysDown_GOBP_TUM], gene_list_LEC_CTRL_TUM, fgseaRes_GOBP_TUM, 
              gseaParam = 0.5)
dev.off()
plotGseaTable(GOBP[topPathwaysDown_GOBP_VEGFC], gene_list_LEC_CTRL_VEGFC, fgseaRes_GOBP_VEGFC, 
              gseaParam = 0.5)
dev.off()
plotGseaTable(KEGG[topPathwaysDown_KEGG_TUM], gene_list_LEC_CTRL_TUM, fgseaRes_KEGG_TUM, 
              gseaParam = 0.5)
dev.off()
plotGseaTable(KEGG[topPathwaysDown_KEGG_VEGFC], gene_list_LEC_CTRL_VEGFC, fgseaRes_KEGG_VEGFC, 
              gseaParam = 0.5)
dev.off()
plotGseaTable(REACTOME[topPathwaysDown_REACTOME_TUM], gene_list_LEC_CTRL_TUM, fgseaRes_REACTOME_TUM, 
              gseaParam = 0.5)
dev.off()
plotGseaTable(REACTOME[topPathwaysDown_REACTOME_VEGFC], gene_list_LEC_CTRL_VEGFC, fgseaRes_REACTOME_VEGFC, 
              gseaParam = 0.5)
dev.off()

# Enrichment plots ----
# Pathways enriched in teLEC
plotEnrichment(GOBP[["GO_RESPONSE_TO_CHEMOKINE"]], gene_list_LEC_CTRL_TUM)
plotEnrichment(GOBP[["GO_STEROL_BIOSYNTHETIC_PROCESS"]], gene_list_LEC_CTRL_TUM)
plotEnrichment(KEGG[["KEGG_CYTOKINE_CYTOKINE_RECEPTOR_INTERACTION"]], gene_list_LEC_CTRL_TUM)
plotEnrichment(KEGG[["KEGG_STEROID_BIOSYNTHESIS"]], gene_list_LEC_CTRL_TUM)
plotEnrichment(REACTOME[["REACTOME_SIGNALING_BY_INTERLEUKINS"]], gene_list_LEC_CTRL_TUM)
plotEnrichment(REACTOME[["REACTOME_CHOLESTEROL_BIOSYNTHESIS"]], gene_list_LEC_CTRL_TUM)
# Pathways enriched in VEGF-LEC
plotEnrichment(GOBP[["GO_CHROMOSOME_SEGREGATION"]], gene_list_LEC_CTRL_VEGFC)
plotEnrichment(GOBP[["GO_CHLORIDE_TRANSPORT"]], gene_list_LEC_CTRL_VEGFC)
plotEnrichment(KEGG[["KEGG_CELL_CYCLE"]], gene_list_LEC_CTRL_VEGFC)
plotEnrichment(KEGG[["KEGG_MAPK_SIGNALING_PATHWAY"]], gene_list_LEC_CTRL_VEGFC)
plotEnrichment(REACTOME[["REACTOME_CELL_CYCLE_CHECKPOINTS"]], gene_list_LEC_CTRL_VEGFC)
plotEnrichment(REACTOME[["REACTOME_ION_CHANNEL_TRANSPORT"]], gene_list_LEC_CTRL_VEGFC)
dev.off()

# Heatmaps of the leading edges (LE) ----
# Represent the log2FC values for the genes that contribute most to the ES

# Pathways enriched in teLEC
# Heatmap for upregulated pathways (positive enrichment score)
LE <- fgseaRes_GOBP_TUM[fgseaRes_GOBP_TUM$pathway == "GO_RESPONSE_TO_CHEMOKINE", 8][[1]][[1]]
# Can be done for other upregulated pathways if needed
# Make a dataframe
df_LE <- as.data.frame(LE)
row.names(df_LE) <- df_LE$LE
head(df_LE)
# Merge with the tables obtained after DE (all genes, logFC, pval, padj)
DE_LEC_CTRL_TUM_all_genes <- read.table("DE_LEC_CTRL_TUM_all_genes.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_TUM_all_genes)
nrow(DE_LEC_CTRL_TUM_all_genes)
DE_LEC_CTRL_VEGFC_all_genes <- read.table("DE_LEC_CTRL_VEGFC_all_genes.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_VEGFC_all_genes)
nrow(DE_LEC_CTRL_VEGFC_all_genes)
merge1 <- merge(df_LE, DE_LEC_CTRL_TUM_all_genes, by = "row.names", all.x = TRUE)
head(merge1)
row.names(merge1) <- merge1$Row.names
merge1 <- merge1[ , c(3,7)]
head(merge1)
colnames(merge1) <- c("logFC_teLEC", "padj_teLEC")
head(merge1)
# Add results obtained for those genes in the VEGF-LEC condition (as comparison)
merge2 <- merge(merge1, DE_LEC_CTRL_VEGFC_all_genes, by = "row.names", all.x = TRUE)
head(merge2)
row.names(merge2) <- merge2$Row.names
merge2 <- merge2[ , c(2:4,8)]
head(merge2)
colnames(merge2) <- c("logFC_teLEC", "padj_teLEC","logFC_VEGFLEC", "padj_VEGFLEC")
head(merge2)
# If padj values > 0.05 (non-significance), set the logFC values to 0 (as these are non-DE genes)
x <- merge2
x$logFC_teLEC <- x$logFC_teLEC * (as.integer(x$padj_teLEC < 0.05))
x$logFC_teLEC
x$logFC_VEGFLEC <- x$logFC_VEGFLEC * (as.integer(x$padj_VEGFLEC < 0.05))
x$logFC_VEGFLEC
# Keep only the log2FC
data_logFC <- x[ , c(1,3)]
head(data_logFC)
colnames(data_logFC) <- c("teLEC", "VEGF-LEC")
head(data_logFC)
# Order by decreasing log2FC in teLEC
data_logFC_ordered <- data_logFC[order(-data_logFC$teLEC), ]
head(data_logFC_ordered)
# Convert to a matrix
data_mat <- data.matrix(data_logFC_ordered)
# Take the first 12 genes of the LE (and transpose the matrix)
LE_top_12_up <- t(head(data_mat, 12))
ncol <- 100
rampbreaks <- seq(-6, 6, length.out = ncol+1)
pheatmap(LE_top_12_up, 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = rampbreaks,
         scale = "none",
         legend = TRUE,
         legend_breaks = c(-6,0,6),
         legend_labels = c(-6,0,6),
         angle_col = 45,
         fontsize_row = 8,
         fontsize_col = 8,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         cellwidth = 30,
         cellheight = 20,
         treeheight_row = 0,
         display_numbers = TRUE,
         na_col = "white")

# Heatmap for downregulated pathways (negative enrichment score)
LE <- fgseaRes_GOBP_TUM[fgseaRes_GOBP_TUM$pathway == "GO_STEROL_BIOSYNTHETIC_PROCESS", 8][[1]][[1]]
# Can be done for other downregulated pathways if needed
# Make a dataframe
df_LE <- as.data.frame(LE)
row.names(df_LE) <- df_LE$LE
head(df_LE)
# Merge with the tables obtained after DE (all genes, logFC, pval, padj)
DE_LEC_CTRL_TUM_all_genes <- read.table("DE_LEC_CTRL_TUM_all_genes.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_TUM_all_genes)
nrow(DE_LEC_CTRL_TUM_all_genes)
DE_LEC_CTRL_VEGFC_all_genes <- read.table("DE_LEC_CTRL_VEGFC_all_genes.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_VEGFC_all_genes)
nrow(DE_LEC_CTRL_VEGFC_all_genes)
merge1 <- merge(df_LE, DE_LEC_CTRL_TUM_all_genes, by = "row.names", all.x = TRUE)
head(merge1)
row.names(merge1) <- merge1$Row.names
merge1 <- merge1[ , c(3,7)]
head(merge1)
colnames(merge1) <- c("logFC_teLEC", "padj_teLEC")
head(merge1)
# Add results obtained for those genes in the VEGF-LEC condition (as comparison)
merge2 <- merge(merge1, DE_LEC_CTRL_VEGFC_all_genes, by = "row.names", all.x = TRUE)
head(merge2)
row.names(merge2) <- merge2$Row.names
merge2 <- merge2[ , c(2:4,8)]
head(merge2)
colnames(merge2) <- c("logFC_teLEC", "padj_teLEC","logFC_VEGFLEC", "padj_VEGFLEC")
head(merge2)
# If padj values > 0.05 (non-significance), set the logFC values to 0 (as these are non-DE genes)
x <- merge2
x$logFC_teLEC <- x$logFC_teLEC * (as.integer(x$padj_teLEC < 0.05))
x$logFC_teLEC
x$logFC_VEGFLEC <- x$logFC_VEGFLEC * (as.integer(x$padj_VEGFLEC < 0.05))
x$logFC_VEGFLEC
# Keep only the log2FC
data_logFC <- x[ , c(1,3)]
head(data_logFC)
colnames(data_logFC) <- c("teLEC", "VEGF-LEC")
head(data_logFC)
# Order by decreasing log2FC in teLEC
data_logFC_ordered <- data_logFC[order(-data_logFC$teLEC), ]
head(data_logFC_ordered)
# Convert to a matrix
data_mat <- data.matrix(data_logFC_ordered)
# Take the last 12 genes of the LE (and transpose the matrix)
LE_top_12_down <- t(tail(data_mat, 12))
ncol <- 100
rampbreaks <- seq(-6, 6, length.out = ncol+1)
pheatmap(LE_top_12_down, 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = rampbreaks,
         scale = "none",
         legend = TRUE,
         legend_breaks = c(-6,0,6),
         legend_labels = c(-6,0,6),
         angle_col = 45,
         fontsize_row = 8,
         fontsize_col = 8,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         cellwidth = 30,
         cellheight = 20,
         treeheight_row = 0,
         display_numbers = TRUE,
         na_col = "white")

dev.off()

# Pathways enriched in VEGF-LEC
# Heatmap for upregulated pathways (positive enrichment score)
LE <- fgseaRes_GOBP_VEGFC[fgseaRes_GOBP_VEGFC$pathway == "GO_CHROMOSOME_SEGREGATION", 8][[1]][[1]]
# Can be done for other upregulated pathways if needed
# Make a dataframe
df_LE <- as.data.frame(LE)
row.names(df_LE) <- df_LE$LE
head(df_LE)
# Merge with the tables obtained after DE (all genes, logFC, pval, padj)
DE_LEC_CTRL_TUM_all_genes <- read.table("DE_LEC_CTRL_TUM_all_genes.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_TUM_all_genes)
nrow(DE_LEC_CTRL_TUM_all_genes)
DE_LEC_CTRL_VEGFC_all_genes <- read.table("DE_LEC_CTRL_VEGFC_all_genes.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_VEGFC_all_genes)
nrow(DE_LEC_CTRL_VEGFC_all_genes)
merge1 <- merge(df_LE, DE_LEC_CTRL_VEGFC_all_genes, by = "row.names", all.x = TRUE)
head(merge1)
row.names(merge1) <- merge1$Row.names
merge1 <- merge1[ , c(3,7)]
head(merge1)
colnames(merge1) <- c("logFC_VEGFLEC", "padj_VEGFLEC")
head(merge1)
# Add results obtained for those genes in the teLEC condition (as comparison)
merge2 <- merge(merge1, DE_LEC_CTRL_TUM_all_genes, by = "row.names", all.x = TRUE)
head(merge2)
row.names(merge2) <- merge2$Row.names
merge2 <- merge2[ , c(2:4,8)]
head(merge2)
colnames(merge2) <- c("logFC_VEGFLEC", "padj_VEGFLEC","logFC_teLEC", "padj_teLEC")
head(merge2)
# If padj values > 0.05 (non-significance), set the logFC values to 0 (as these are non-DE genes)
x <- merge2
x$logFC_teLEC <- x$logFC_teLEC * (as.integer(x$padj_teLEC < 0.05))
x$logFC_teLEC
x$logFC_VEGFLEC <- x$logFC_VEGFLEC * (as.integer(x$padj_VEGFLEC < 0.05))
x$logFC_VEGFLEC
# Keep only the log2FC
data_logFC <- x[ , c(1,3)]
head(data_logFC)
colnames(data_logFC) <- c("VEGF-LEC", "teLEC")
head(data_logFC)
# Order by decreasing log2FC in VEGF-LEC
data_logFC_ordered <- data_logFC[order(-data_logFC$`VEGF-LEC`), ]
head(data_logFC_ordered)
# Convert to a matrix
data_mat <- data.matrix(data_logFC_ordered)
# Take the first 12 genes of the LE (and transpose the matrix)
LE_top_12_up <- t(head(data_mat, 12))
ncol <- 100
rampbreaks <- seq(-6, 6, length.out = ncol+1)
pheatmap(LE_top_12_up, 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = rampbreaks,
         scale = "none",
         legend = TRUE,
         legend_breaks = c(-6,0,6),
         legend_labels = c(-6,0,6),
         angle_col = 45,
         fontsize_row = 8,
         fontsize_col = 8,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         cellwidth = 30,
         cellheight = 20,
         treeheight_row = 0,
         display_numbers = TRUE,
         na_col = "white")

# Heatmap for downregulated pathways (negative enrichment score)
LE <- fgseaRes_GOBP_VEGFC[fgseaRes_GOBP_VEGFC$pathway == "GO_CHLORIDE_TRANSPORT", 8][[1]][[1]]
# Can be done for other downregulated pathways if needed
# Make a dataframe
df_LE <- as.data.frame(LE)
row.names(df_LE) <- df_LE$LE
head(df_LE)
# Merge with the tables obtained after DE (all genes, logFC, pval, padj)
DE_LEC_CTRL_TUM_all_genes <- read.table("DE_LEC_CTRL_TUM_all_genes.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_TUM_all_genes)
nrow(DE_LEC_CTRL_TUM_all_genes)
DE_LEC_CTRL_VEGFC_all_genes <- read.table("DE_LEC_CTRL_VEGFC_all_genes.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_VEGFC_all_genes)
nrow(DE_LEC_CTRL_VEGFC_all_genes)
merge1 <- merge(df_LE, DE_LEC_CTRL_VEGFC_all_genes, by = "row.names", all.x = TRUE)
head(merge1)
row.names(merge1) <- merge1$Row.names
merge1 <- merge1[ , c(3,7)]
head(merge1)
colnames(merge1) <- c("logFC_VEGFLEC", "padj_VEGFLEC")
head(merge1)
# Add results obtained for those genes in the teLEC condition (as comparison)
merge2 <- merge(merge1, DE_LEC_CTRL_TUM_all_genes, by = "row.names", all.x = TRUE)
head(merge2)
row.names(merge2) <- merge2$Row.names
merge2 <- merge2[ , c(2:4,8)]
head(merge2)
colnames(merge2) <- c("logFC_VEGFLEC", "padj_VEGFLEC","logFC_teLEC", "padj_teLEC")
head(merge2)
# If padj values > 0.05 (non-significance), set the logFC values to 0 (as these are non-DE genes)
x <- merge2
x$logFC_teLEC <- x$logFC_teLEC * (as.integer(x$padj_teLEC < 0.05))
x$logFC_teLEC
x$logFC_VEGFLEC <- x$logFC_VEGFLEC * (as.integer(x$padj_VEGFLEC < 0.05))
x$logFC_VEGFLEC
# Keep only the log2FC
data_logFC <- x[ , c(1,3)]
head(data_logFC)
colnames(data_logFC) <- c("VEGF-LEC", "teLEC")
head(data_logFC)
# Order by decreasing log2FC in VEGF-LEC
data_logFC_ordered <- data_logFC[order(-data_logFC$`VEGF-LEC`), ]
head(data_logFC_ordered)
# Convert to a matrix
data_mat <- data.matrix(data_logFC_ordered)
# Take the last 12 genes of the LE (and transpose the matrix)
LE_top_12_down <- t(tail(data_mat, 12))
ncol <- 100
rampbreaks <- seq(-6, 6, length.out = ncol+1)
pheatmap(LE_top_12_down, 
         color = colorRampPalette(c("blue", "white", "red"))(100),
         breaks = rampbreaks,
         scale = "none",
         legend = TRUE,
         legend_breaks = c(-6,0,6),
         legend_labels = c(-6,0,6),
         angle_col = 45,
         fontsize_row = 8,
         fontsize_col = 8,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         cellwidth = 30,
         cellheight = 20,
         treeheight_row = 0,
         display_numbers = TRUE,
         na_col = "white")

dev.off()

# KEGG over-representation test  ----

# Create gene lists
DE_LEC_CTRL_TUM_all_genes <- read.table("DE_LEC_CTRL_TUM_all_genes.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_TUM_all_genes)
DE_LEC_CTRL_VEGFC_all_genes <- read.table("DE_LEC_CTRL_VEGFC_all_genes.txt", header=TRUE, sep="\t")
head(DE_LEC_CTRL_VEGFC_all_genes)
gene_list_LEC_CTRL_TUM = DE_LEC_CTRL_TUM_all_genes$logFC
gene_list_LEC_CTRL_VEGFC = DE_LEC_CTRL_VEGFC_all_genes$logFC
names(gene_list_LEC_CTRL_TUM) = row.names(DE_LEC_CTRL_TUM_all_genes)
names(gene_list_LEC_CTRL_VEGFC) = row.names(DE_LEC_CTRL_VEGFC_all_genes)
str(gene_list_LEC_CTRL_TUM)
str(gene_list_LEC_CTRL_VEGFC)
gene_list_LEC_CTRL_TUM = sort(gene_list_LEC_CTRL_TUM, decreasing = TRUE)
gene_list_LEC_CTRL_VEGFC = sort(gene_list_LEC_CTRL_VEGFC, decreasing = TRUE)
gene_list_LEC_CTRL_TUM = gene_list_LEC_CTRL_TUM[!duplicated(names(gene_list_LEC_CTRL_TUM))]
gene_list_LEC_CTRL_VEGFC = gene_list_LEC_CTRL_VEGFC[!duplicated(names(gene_list_LEC_CTRL_VEGFC))]
str(gene_list_LEC_CTRL_TUM)
str(gene_list_LEC_CTRL_VEGFC)
head(gene_list_LEC_CTRL_TUM)
head(gene_list_LEC_CTRL_VEGFC)
gene_IDs_LEC_CTRL_TUM <- names(gene_list_LEC_CTRL_TUM)
gene_IDs_LEC_CTRL_VEGFC <- names(gene_list_LEC_CTRL_VEGFC)
head(gene_IDs_LEC_CTRL_TUM)
head(gene_IDs_LEC_CTRL_VEGFC)

# We need entrez IDs for this KEGG over-representation test
# Find the corresponding entrezIDs
# Using biomaRt package
listMarts()
ensembl <- useMart("ensembl")
datasets <- listDatasets(ensembl)
head(datasets)
ensembl = useDataset("hsapiens_gene_ensembl",mart=ensembl)

entrezHUGOGid_LEC_CTRL_TUM <- getBM(attributes = c('entrezgene_id','hgnc_symbol'), 
                                    filters = 'hgnc_symbol', 
                                    values = gene_IDs_LEC_CTRL_TUM, 
                                    mart = ensembl)
entrezHUGOGid_LEC_CTRL_VEGFC <- getBM(attributes = c('entrezgene_id','hgnc_symbol'), 
                                      filters = 'hgnc_symbol', 
                                      values = gene_IDs_LEC_CTRL_VEGFC, 
                                      mart = ensembl)

head(entrezHUGOGid_LEC_CTRL_TUM)
head(entrezHUGOGid_LEC_CTRL_VEGFC)
colnames(entrezHUGOGid_LEC_CTRL_TUM) <- c("entrez_id", "geneID")
colnames(entrezHUGOGid_LEC_CTRL_VEGFC) <- c("entrez_id", "geneID")
head(entrezHUGOGid_LEC_CTRL_TUM)
head(entrezHUGOGid_LEC_CTRL_VEGFC)

# Add a "geneID" column to tables of DE genes (prepare for merging datasets)
head(DE_LEC_CTRL_TUM_all_genes)
head(DE_LEC_CTRL_VEGFC_all_genes)
DE_LEC_CTRL_TUM_all_genes$geneID <- rownames(DE_LEC_CTRL_TUM_all_genes)
DE_LEC_CTRL_VEGFC_all_genes$geneID <- rownames(DE_LEC_CTRL_VEGFC_all_genes)
head(DE_LEC_CTRL_TUM_all_genes)
head(DE_LEC_CTRL_VEGFC_all_genes)

# Merge entrezIDs and tables of DE genes
DE_LEC_CTRL_TUM_with_entrezID <- merge(entrezHUGOGid_LEC_CTRL_TUM, DE_LEC_CTRL_TUM_all_genes, by = "geneID")
DE_LEC_CTRL_VEGFC_with_entrezID <- merge(entrezHUGOGid_LEC_CTRL_VEGFC, DE_LEC_CTRL_VEGFC_all_genes, by = "geneID")
head(DE_LEC_CTRL_TUM_with_entrezID)
head(DE_LEC_CTRL_VEGFC_with_entrezID)
# Remove NA
DE_LEC_CTRL_TUM_with_entrezID <- na.exclude(DE_LEC_CTRL_TUM_with_entrezID)
DE_LEC_CTRL_VEGFC_with_entrezID <- na.exclude(DE_LEC_CTRL_VEGFC_with_entrezID)
head(DE_LEC_CTRL_TUM_with_entrezID)
head(DE_LEC_CTRL_VEGFC_with_entrezID)
# Order by decreasing log2FC values
DE_LEC_CTRL_TUM_with_entrezID_by_logFC <- DE_LEC_CTRL_TUM_with_entrezID[order(-DE_LEC_CTRL_TUM_with_entrezID$logFC), ]
DE_LEC_CTRL_VEGFC_with_entrezID_by_logFC <- DE_LEC_CTRL_VEGFC_with_entrezID[order(-DE_LEC_CTRL_VEGFC_with_entrezID$logFC), ]
head(DE_LEC_CTRL_TUM_with_entrezID_by_logFC)
head(DE_LEC_CTRL_VEGFC_with_entrezID_by_logFC)

# Select significant DE genes
sigGenes_LEC_CTRL_TUM <- DE_LEC_CTRL_TUM_with_entrezID_by_logFC$entrez_id[DE_LEC_CTRL_TUM_with_entrezID_by_logFC$padj < 0.05 & 
                                                                            abs(DE_LEC_CTRL_TUM_with_entrezID_by_logFC$logFC) > 1 ]
sigGenes_LEC_CTRL_VEGFC <- DE_LEC_CTRL_VEGFC_with_entrezID_by_logFC$entrez_id[DE_LEC_CTRL_VEGFC_with_entrezID_by_logFC$padj < 0.05 & 
                                                                                abs(DE_LEC_CTRL_VEGFC_with_entrezID_by_logFC$logFC) > 1 ]

sigGenes_LEC_CTRL_TUM <- na.exclude(sigGenes_LEC_CTRL_TUM)
sigGenes_LEC_CTRL_VEGFC <- na.exclude(sigGenes_LEC_CTRL_VEGFC)
head(sigGenes_LEC_CTRL_TUM)
head(sigGenes_LEC_CTRL_VEGFC)

# KEGG enrichment analysis
KEGG_enrichment_LEC_CTRL_TUM <- enrichKEGG(gene = sigGenes_LEC_CTRL_TUM, organism = 'hsa')
KEGG_enrichment_LEC_CTRL_TUM[, c(2:6) ]
KEGG_enrichment_LEC_CTRL_VEGFC <- enrichKEGG(gene = sigGenes_LEC_CTRL_VEGFC, organism = 'hsa')
KEGG_enrichment_LEC_CTRL_VEGFC[, c(2:6) ]

