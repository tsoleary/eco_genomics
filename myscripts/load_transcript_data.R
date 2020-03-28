# Loading Red Spruce transcriptomics data --------------------------------------

# Packages
require(tidyverse)
require(DESeq2)

# Load counts data -------------------------------------------------------------

# Set working directory
setwd(here::here("myresults/transcript_counts"))

# Import and round data from Salmon
counts <- round(read.table("RS_cds2kb_countsMatrix.txt", header = TRUE))

# Load the metadata
metadata <- read_delim("RS_samples.txt", delim = "\t",
                       col_types = "cffff")

# Filter to only have the Day 10 samples 
metadata <- metadata %>%
  filter(day == "10")

counts <- counts %>%
  rownames_to_column("gene") %>%
  select(c("gene", metadata$sample)) %>%
  column_to_rownames("gene")


# Run DESeq2 -------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = 
                                ~ climate + treatment + climate:treatment)

# Set the group to be compared against
dds$group <- relevel(dds$climate, ref = "CW")

# Filter out genes with few reads
dds <- dds[rowSums(counts(dds)) > 30]

# Run DESeq
dds <- DESeq(dds)

# List the results you've generated
resultsNames(dds)

res_climate <- results(dds, contrast = c("climate", "HD", "CW"), alpha = 0.05)
res_D_C <- results(dds, contrast = c("treatment", "D", "C"), alpha = 0.05)
res_H_C <- results(dds, contrast = c("treatment", "H", "C"), alpha = 0.05)
#res_int <- results(dds, name = "climateCW.treatmentH")
#res_int <- results(dds, name = "climateCW.treatmentD")

# summary(res_climate)
# summary(res_D_C)
# summary(res_H_C)

# Transform for PCA ------------------------------------------------------------
vsd <- vst(dds)
pc_data <- plotPCA(vsd, 
                   ntop = 100000,
                   intgroup = c("climate", "treatment"), 
                   returnData = TRUE)
percentVar <- round(100 * attr(pc_data, "percentVar"))

# reordering the factors
pc_data$treatment <- factor(pc_data$treatment, 
                         levels = c("C", "H", "D"), 
                         labels = c("C", "H", "D"))

# Combining factors ------------------------------------------------------------
metadata <- metadata %>%
  mutate(group = paste(climate, treatment, sep = "_"))

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ group)

# Set the group to be compared against CW comparisons -----
dds$group <- relevel(dds$group, ref = "CW_C")

# Filter out genes with few reads
dds <- dds[rowSums(counts(dds)) > 30]

# Run DESeq
dds <- DESeq(dds)

# List the results you've generated
resultsNames(dds)

res_CW_D <- results(dds, contrast = c("group", "CW_D", "CW_C"), alpha = 0.05)
res_CW_H <- results(dds, contrast = c("group", "CW_H", "CW_C"), alpha = 0.05)

# summary(res_CW_D)
# summary(res_CW_H)

# Set the group to be compared against HD comparisons -----
dds$group <- relevel(dds$group, ref = "HD_C")

# Filter out genes with few reads
dds <- dds[rowSums(counts(dds)) > 30]

# Run DESeq
dds <- DESeq(dds)

# List the results you've generated
resultsNames(dds)

res_HD_D <- results(dds, contrast = c("group", "HD_D", "HD_C"), alpha = 0.05)
res_HD_H <- results(dds, contrast = c("group", "HD_H", "HD_C"), alpha = 0.05)

# summary(res_HD_D)
# summary(res_HD_H)

