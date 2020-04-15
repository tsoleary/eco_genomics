# ------------------------------------------------------------------------------
# P. rubens transcriptomics for GWAS-Env project
# April 15, 2020
# TS O'Leary
# ------------------------------------------------------------------------------

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


# Run DESeq2 -------------------------------------------------------------------
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = 
                                ~ climate + treatment + day + 
                                climate:treatment + treatment:day + climate:day)

# # Set the group to be compared against
dds$group <- relevel(dds$day, ref = "0")

# Filter out genes with few reads
dds <- dds[rowSums(counts(dds)) > 76]

# Run DESeq
dds <- DESeq(dds)
setwd(here::here("myresults/gwas_env"))
saveRDS(dds, "p_rubens_clim_treat_day_all_dds.rds")

# List the results you've generated
resultsNames(dds)

res_climate <- results(dds, contrast = c("climate", "HD", "CW"), alpha = 0.05)
res_D_C <- results(dds, contrast = c("treatment", "D", "C"), alpha = 0.05)
res_H_C <- results(dds, contrast = c("treatment", "H", "C"), alpha = 0.05)
res_5_0 <- results(dds, contrast = c("day", "5", "0"), alpha = 0.05)
res_10_0 <- results(dds, contrast = c("day", "10", "0"), alpha = 0.05)

# Summary stats
summary(res_climate)
summary(res_D_C)
summary(res_H_C)
summary(res_5_0)
summary(res_10_0)


setwd(here::here("myresults/gwas_env"))
saveRDS(res_climate, "p_rubens_res_climate.rds")
saveRDS(res_D_C, "p_rubens_res_D_C.rds")
saveRDS(res_H_C, "p_rubens_res_H_C.rds")
saveRDS(res_5_0, "p_rubens_res_5_0.rds")
saveRDS(res_10_0, "p_rubens_res_10_0.rds")

# Combining factors to run pairwise tests --------------------------------------
metadata <- metadata %>%
  mutate(group = paste(climate, treatment, day, sep = "_"))


# Filter to include only samples with three or more replicates
met_filt <- metadata %>%
  dplyr::group_by(group) %>%
  dplyr::count() %>%
  dplyr::filter(n >= 3)

metadata <- metadata %>%
  filter(group %in% met_filt$group)

counts <- counts %>%
  rownames_to_column("gene") %>%
  select(c("gene", metadata$sample)) %>%
  column_to_rownames("gene")


# Set up dds object with counts, metadata, and experimental design
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = metadata,
                              design = ~ group)

# Filter out genes with few reads
dds <- dds[rowSums(counts(dds)) > 73]


# Cool & Wet: Control time series ie CW_C_0 comparisons ------------------------
dds$group <- relevel(dds$group, ref = "CW_C_0")

# Run DESeq
dds <- DESeq(dds)

# List the results you've generated
resultsNames(dds)

res_CW_C_5_0 <- results(dds, 
                        contrast = c("group", "CW_C_5", "CW_C_0"), 
                        alpha = 0.05)
res_CW_C_10_0 <- results(dds, 
                         contrast = c("group", "CW_C_10", "CW_C_0"), 
                         alpha = 0.05)

summary(res_CW_C_5_0)
summary(res_CW_C_10_0)

setwd(here::here("myresults/gwas_env"))
saveRDS(dds, "p_rubens_CW_C_0_dds.rds")
saveRDS(res_CW_C_5_0, "p_rubens_res_CW_C_5_0.rds")
saveRDS(res_CW_C_10_0, "p_rubens_res_CW_C_10_0.rds")


# Cool & Wet: Heat only stress time series ie CW_H_0 comparisons ---------------
dds$group <- relevel(dds$group, ref = "CW_H_0")

# Run DESeq
dds <- DESeq(dds)

# List the results you've generated
resultsNames(dds)

res_CW_H_5_0 <- results(dds, 
                        contrast = c("group", "CW_H_5", "CW_H_0"), 
                        alpha = 0.05)
res_CW_H_10_0 <- results(dds, 
                         contrast = c("group", "CW_H_10", "CW_H_0"), 
                         alpha = 0.05)


summary(res_CW_H_5_0)
summary(res_CW_H_10_0)


setwd(here::here("myresults/gwas_env"))
saveRDS(dds, "p_rubens_CW_H_0_dds.rds")
saveRDS(res_CW_H_5_0, "p_rubens_res_CW_H_5_0.rds")
saveRDS(res_CW_H_10_0, "p_rubens_res_CW_H_10_0.rds")


# Cool & Wet: Heat & Drought stress time series ie CW_D_0 comparisons ----------
dds$group <- relevel(dds$group, ref = "CW_D_0")

# Run DESeq
dds <- DESeq(dds)

# List the results you've generated
resultsNames(dds)

# CW_D_5 did not have three or more reps
# res_CW_D_5_0 <- results(dds, 
#                         contrast = c("group", "CW_D_5", "CW_D_0"), 
#                         alpha = 0.05)
res_CW_D_10_0 <- results(dds, 
                         contrast = c("group", "CW_D_10", "CW_D_0"), 
                         alpha = 0.05)


# summary(res_CW_D_5_0)
summary(res_CW_D_10_0)


setwd(here::here("myresults/gwas_env"))
saveRDS(dds, "p_rubens_CW_D_0_dds.rds")
# saveRDS(res_CW_D_5_0, "p_rubens_res_CW_D_5_0.rds")
saveRDS(res_CW_D_10_0, "p_rubens_res_CW_D_10_0.rds")

