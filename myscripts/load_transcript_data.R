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
                       col_types = cols(
                         sample = col_character(),
                         pop = col_factor(),
                         treatment = col_factor(),
                         day = col_factor(),
                         climate = col_factor()))

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

summary(res_climate)
summary(res_D_C)
summary(res_H_C)

# splitting the groups up with paste and looking at it that way?....
metadata <- metadata %>%
  mutate(group = paste(climate, treatment, sep = "_"))

