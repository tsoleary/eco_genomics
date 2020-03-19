# Transcriptomics in R

# Set your working directory
setwd(here::here("myresults/transcript_counts"))

# Packages 
library(DESeq2)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)
library(ggpubr)
library(wesanderson)
library(vsn)
library(pheatmap)

# Import data ------------------------------------------------------------------
countsTable <- read.table("RS_cds2kb_countsMatrix.txt", 
                          header=TRUE, 
                          row.names=1)
head(countsTable)
dim(countsTable)
# Need to round because DESeq wants only integers
countsTableRound <- round(countsTable)
head(countsTableRound)

# Import the samples description table 
# links each sample to factors of the experimental design
# Need the colClasses otherwise imports "day" as numeric which DESeq doesn't do
conds <- read.delim("RS_samples.txt", 
                    header = TRUE, 
                    stringsAsFactors = TRUE, 
                    row.names = 1, 
                    colClasses = c('factor', 'factor', 'factor', 'factor'))
head(conds)
dim(conds)


# Read summary statistics and visualization ------------------------------------
colSums(countsTableRound)
mean(colSums(countsTableRound))
# bar plot to vizualize the sum of reads
barplot(colSums(countsTableRound), 
        las = 3, 
        cex.names = 0.5,
        names.arg = substring(colnames(countsTableRound), 1, 13))
abline(h = mean(colSums(countsTableRound)), col = "blue", lwd = 2)

# Avergave number of counts per gene
rowSums(countsTableRound)
apply(countsTableRound, 2, mean)

# Create a DESeq object --------------------------------------------------------
# define the experimental design here with the tilde

dds <- DESeqDataSetFromMatrix(countData = countsTableRound,
                              colData = conds,
                              design = ~ climate + day + treatment)

# Filter out genes with few reads
# 76 was chosen to create an average of 1 read per sample
dds <- dds[rowSums(counts(dds)) > 76]

# Run the DESeq model to test for differential gene expression -----------------
# DESeq does these three things all in this function
# 1) estimate size factors (per sample) 
# 2) estimate dispersion (per gene)
# 3) run negative binomial glm
dds <- DESeq(dds)

# List the results you've generated
resultsNames(dds)

# Order and list and summarize results from specific contrasts -----------------
# Set your adjusted p-value cutoff
# Can make summary tables of the number of genes differentially expressed 
# up/down-regulated) for each contrast

res_treatCD <- results(dds, contrast = c("treatment", "C", "D"), alpha = 0.05)
res_treatCD <- res_treatCD[order(res_treatCD$padj), ]
summary(res_treatCD)

# Data visualization -----------------------------------------------------------
# MA plot
ggmaplot(res_treatCD, 
         main = expression("Hot & Dry" %->% "Control"),
         fdr = 0.01, fc = 0, size = 1,
         palette = c("#B31B21", "#1465AC", "#A6A6A680"),
         legend = "top", 
         top = 0,
         # genenames = as.vector(res_treatCD$gene),
         # select.top.method = "padj",
         # font.label = c("bold", 11), 
         # label.rectangle = TRUE,
         font.legend = "bold",
         font.main = "bold",
         ggtheme = ggplot2::theme_classic()) + 
  theme(legend.text = element_text(size = 12))

# PCA
vsd <- vst(dds, blind = FALSE, nsub = 10000)
data <- plotPCA(vsd, 
                ntop = 10000,
                intgroup = c("climate", "treatment", "day"), 
                returnData = TRUE)
percentVar <- round(100 * attr(data, "percentVar"))

# reordering the factors
data$treatment <- factor(data$treatment, 
                         levels = c("C", "H", "D"), 
                         labels = c("C", "H", "D"))
data$day <- factor(data$day, 
                   levels = c("0", "5", "10"), 
                   labels = c("0","5","10"))

ggplot(data, aes(PC1, PC2, color = climate, shape=treatment)) +
  geom_point(size = 4, alpha = 0.85) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  theme_minimal()


ggplot(data) +
  geom_histogram(aes(x = PC1, fill = day), 
                 bins = 20, color = "grey50", position = "dodge") +
  theme_classic()

# Counts of specific top gene! 
# important validatition that the normalization, model is working
d <- plotCounts(dds, 
                gene = "MA_10426407g0030", 
                intgroup = (c("treatment","climate", "day")), 
                returnData = TRUE)

ggplot(d, aes(x = treatment, 
                   y = count, 
                   shape = climate, 
                   color = day)) + 
  geom_jitter(width = 0.3, size = 3) +
  scale_x_discrete(limits = c("C", "H", "D")) +
  theme_classic() + 
  theme(text = element_text(size = 20), 
        panel.grid.major = element_line(colour = "grey"))

# Heatmap of top 20 genes sorted by pvalue
topgenes <- head(rownames(res_treatCD), 20)
mat <- assay(vsd)[topgenes, ]
mat <- mat - rowMeans(mat)
df <- as.data.frame(colData(dds)[, c("treatment", "climate", "day")])
pheatmap(mat, 
         annotation_col = df, 
         show_colnames = FALSE)



