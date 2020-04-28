# ------------------------------------------------------------------------------
# Weighted gene correlation network analysis on P.rubens transcriptomics data
# April 22, 2020
# Katelynn Warner, Ben Camber, & TS O'Leary
# ------------------------------------------------------------------------------

# Load packages 
require(WGCNA)
require(DESeq2)
require(tidyverse)

# Source rds objects into the environment
# In this case, I think using the dds with the full model design is best
setwd(here::here("myresults/gwas_env"))
dds <- readRDS("p_rubens_clim_treat_day_all_dds.rds")


# Extract the normalized counts from the dds object for WGCNA analysis
norm_counts <- counts(dds, normalized = TRUE)

# WGCNA wants this normalized counts matrix transposed
t_norm_counts <- as.data.frame(t(norm_counts))
names(t_norm_counts) <- rownames(norm_counts)

# Check for genes and samples with too many missing values (OK if TRUE)
goodSamplesGenes(t_norm_counts, verbose = 3)$allOK


# Cluster the samples together to determine outliers ---------------------------
# Plot the clusters and save to a pdf
sampleTree <- hclust(dist(t_norm_counts), method = "average")
setwd(here::here("myresults/gwas_env"))
pdf(file = "sampleClustering.pdf", width = 12, height = 9)
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, 
     main = "Sample clustering to detect outliers", 
     sub = "", 
     xlab = "", 
     cex.lab = 1.5,
     cex.axis = 1.5, 
     cex.main = 2)
abline(h = 5000000, col = "red");
dev.off()


# Cut the tree to remove the outlier samples from the dataset
clust <- cutreeStatic(sampleTree, cutHeight = 5000000, minSize = 10)
# clust == 1 corresponds to which samples are to be kept
t_norm_counts <- t_norm_counts[clust == 1, ]
nGenes <- ncol(t_norm_counts)
nSamples <- nrow(t_norm_counts)

# Read in metadata & convert to factors to numeric equivalents & rm pop col
setwd(here::here("myresults/transcript_counts"))
metadata <- read.table("RS_samples.txt", header = TRUE) %>%
  mutate(treatment = case_when(treatment == "C" ~ 1, 
                               treatment == "H" ~ 2,
                               treatment == "D" ~ 3),
         climate = case_when(climate == "CW" ~ 1,
                             climate == "HD" ~ 2)) %>%
  dplyr::select(-pop)

sampleTreeCut <- hclust(dist(t_norm_counts), 
                        method = "average") 

# Filter metadata to have only the samples kept after cutting the tree
metadata_filt <- metadata %>% 
  filter(sample %in% rownames(t_norm_counts)) %>%
  column_to_rownames("sample")

# Plot the cluster dendrogram with the factors labeled with colors
traitColors <- numbers2colors(metadata_filt, 
                              signed = FALSE)
setwd(here::here("myresults/gwas_env"))
pdf(file = "sampleClusteringCutWithFactors.pdf", width = 12, height = 9)
plotDendroAndColors(sampleTreeCut, traitColors,
                    groupLabels = metadata_filt,
                    main = "Sample dendrogram and trait heatmap")
dev.off()


# Network Contruction and module creation --------------------------------------

# Choose a set of soft-thresholding powers
powers <- c(1:10, seq(from = 12, to = 20, by = 2))

# Run the network 
sft <- pickSoftThreshold(t_norm_counts, 
                         powerVector = powers, 
                         verbose = 5)

saveRDS("p_rubens_expr_net_sft.rds")

# Plots
sizeGrWindow(9, 5)
par(mfrow = c(1, 2))
cex1 = 0.9

# Scale-free topology fit indes as a function of the soft-thresholding power
plot(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2], 
     xlab = "Soft Threshold (power)", 
     ylab =  "Scale Free Topology Model Fit, signed R^2", 
     type = "n", 
     main = paste("Scale independence"))

text(sft$fitIndices[, 1], 
     -sign(sft$fitIndices[, 3])*sft$fitIndices[, 2], 
     labels = powers, 
     cex = 1, 
     col = "red")

abline(h = 0.9, col = "red")

# mean connectivity as a function of the soft-threstholding power

plot(sft$fitIndices[, 1], 
     sft$fitIndices[, 5], 
     xlab = "Soft Threshold (power)", 
     ylab = "Mean Connectivity", 
     type = "n", 
     main = paste("Mean connectivity"))
text(sft$fitIndices[, 1], 
     sft$fitIndices[, 5], 
     labels = powers, 
     cex = cex1, 
     col = "red")

# Power 6 is the best description for our model --------------------------------
# This net object created from this code is saved in the file created below
# cor <- WGCNA::cor
# net <- blockwiseModules(t_norm_counts, 
#                         power = 6, TOMType = "unsigned", 
#                         minModuleSize = 30, reassignThreshold = 0, 
#                         mergeCutHeight = 0.25, numericLabels = TRUE, 
#                         pamRespectsDendro = FALSE, saveTOMs = TRUE, 
#                         saveTOMFileBase = "p_rubens_TOM", verbose=3)
# setwd(here::here("myresults/gwas_env"))
# saveRDS(net, "p_rubens_wgcna_net.rds")

# Read in the net object from its .rds
net <- readRDS(here::here("myresults/gwas_env/p_rubens_wgcna_net.rds"))

# Plot the dendrogram with the module colors underneath ------------------------
mergedColors <- labels2colors(net$colors)
plotDendroAndColors(net$dendrograms[[1]],
                    mergedColors[net$blockGenes[[1]]],
                    "Module colors", 
                    dendroLabels = FALSE, 
                    hang = 0.03, 
                    addGuide = TRUE, 
                    guideHang = 0.05)

moduleLabels <- net$colors
moduleColors <- labels2colors(net$colors)
MEs <- net$MEs
geneTree <- net$dendrograms[[1]]

# Define numbers of genes and samples ------------------------------------------

nGenes <- ncol(t_norm_counts)
nSamples <- nrow(t_norm_counts)

MEs0 <- moduleEigengenes(t_norm_counts, moduleColors)$eigengenes
MEs <- orderMEs(MEs0)
moduleTraitCor <- cor(MEs, metadata_filt, use = "p")

moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)


textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", 
                    signif(moduleTraitPvalue,1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)

n_genes_mod <- tibble("names" = names(gene_mod_list),
                      "n_genes" = lengths(gene_mod_list))

lab_df <- tibble("names" = names(MEs)) %>%
  mutate(names = str_replace_all(names, "ME", "")) %>%
  full_join(n_genes_mod) %>%
  mutate(lab = paste0(names, " (", n_genes, ")"))


  


par(mar = c(5.1, 10, 4.1, 2.1))
labeledHeatmap(Matrix = moduleTraitCor, 
               xLabels= names(metadata_filt),
               yLabels = names(MEs), 
               ySymbols = lab_df$lab, 
               colorLabels = FALSE, 
               colors = blueWhiteRed(50), 
               textMatrix = textMatrix, 
               setStdMargins = FALSE, 
               cex.text = 0.5, 
               zlim = c(-1, 1), 
               main = paste("Module-trait relationships"))

# I haven't really looked at or cleaned this shit...
# module <- "MEdarkred"
# moduleGenes <- moduleColors == module
# geneModuleMembership <- as.data.frame(cor(t_norm_counts, MEs, use = "p"))
# MMPvalue <- as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))

# names(geneModuleMembership) <- paste("MM", modNames, sep="")
# names(MMPvalue) = paste("p.MM", modNames, sep="");
# geneTraitSignificance <- as.data.frame(cor(t_norm_counts, day, use = "p"))
# GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
# names(geneTraitSignificance) = paste("GS.", names(day), sep="")
# names(GSPvalue) = paste("p.GS.", names(day), sep="")
# 
# module <- "plum1"
# column <- match(module, modNames);
# moduleGenes <- moduleColors == module;
# sizeGrWindow(7, 7);
# par(mfrow = c(1,1));
# verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
#                    abs(geneTraitSignificance[moduleGenes, 1]),
#                    xlab = paste("Module Membership in", module, "module"),
#                    ylab = "Gene significance for day",
#                    main = paste("Module membership vs. gene significance\n"),
#                    cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "black")
# 
# geneModuleMembership$MMplum1      

# Get vector of genes belonging to each module ---------------------------------
gene_mod_list <- vector(mode = "list", 
                        length = length(unique(moduleColors)))
names(gene_mod_list) <- unique(moduleColors)

saveRDS(gene_mod_list, "gene_mod_list.rds")

for (i in seq_along(unique(moduleColors))) {
  gene_mod_list[[unique(moduleColors)[i]]] <- 
    names(t_norm_counts)[moduleColors == unique(moduleColors)[i]]
}

# Filter only the pvals < 0.01 to get the sig modules --------------------------
sig_mod <- as_tibble(moduleTraitPvalue, rownames = "mod") %>%
  pivot_longer(-mod, names_to = "factor", values_to = "p_val") %>%
  filter(p_val < 0.01)

# Make a data.frame with module, gene, and significance values -----------------
day <- as.data.frame(metadata_filt$day)
names(day) <- "day"
treatment <- as.data.frame(metadata_filt$treatment)
names(treatment) <- "treatment"


geneTraitSignificance_day <- as.data.frame(cor(t_norm_counts, 
                                               day, 
                                               use = "p"))

geneTraitSignificance_treat <- as.data.frame(cor(t_norm_counts, 
                                                 treatment, 
                                                 use = "p"))

sig_genes <- as_tibble(list("mod" = moduleColors, 
                            "genes" = names(t_norm_counts),
                            "sig_day" = geneTraitSignificance_day$day,
                            "sig_treat" = geneTraitSignificance_treat$treat))


# Filter sig_genes only in magenta -----
sig_genes %>%
  filter(mod %in% 
           unique(str_replace(sig_mod$mod, "ME", ""))) %>%
  arrange(mod, desc(sig_day)) %>%
  filter(mod == "darkgrey")

# See what the counts patterns look like in DESeq ------------------------------

vsd <- vst(dds)
df <- as.data.frame(colData(dds)[, c("treatment", "day")])

keepSamples <- rownames(metadata_filt)

# Heatmap
sig_mod_colors <- c("darkgrey", "magenta", "darkred", "red")
setwd(here::here("myresults/gwas_env"))

for (color in sig_mod_colors) {
  topgenes <- gene_mod_list[[color]]
  mat <- assay(vsd)[topgenes, keepSamples]
  mat <- mat - rowMeans(mat)
  pdf(paste0(color, "_heatmap.pdf"), 
      width = 5, height = 10)
  pheatmap::pheatmap(mat, 
                     annotation_col = df, 
                     show_rownames = FALSE,
                     show_colnames = FALSE)
  dev.off()

}



# plot all genes ---------------------------------------------------------------
setwd(here::here("myresults/gwas_env"))
dds <- readRDS("p_rubens_clim_treat_day_all_dds.rds")

# ------------------------------------------------------------------------------
# Function: plot_gene_counts_day
# Description: Plots the normalized counts of a specific gene 
# Inputs: dds DESeq2 object and gene character string
# Outputs: ggplot

plot_gene_counts_day <- function(dds, gene) {
  
  ###subset(dds, rownames(dds@colData) %in% rownames(metadata_filt))
  
  d <- plotCounts(dds, 
                  gene = gene, 
                  intgroup = (c("treatment","climate", "day")), 
                  returnData = TRUE)
  
  ggplot(d, aes(x = day, 
                y = count, 
                shape = climate, 
                color = treatment)) + 
    geom_jitter(width = 0.3, size = 3, alpha = 0.8) +
    scale_x_discrete(limits = c("0", "5", "10")) +
    labs(title = gene,
         y = "Normalized Counts",
         x = "Day") +
    theme_classic() + 
    scale_color_manual(name = "Treatment",
                       values = c("grey50", "blue", "firebrick"),
                       limits = c("C", "H", "D"),
                       labels = c("Control", "Heat", "Heat & Drought")) + 
    scale_shape_manual(name = "Climate",
                       values = c(19, 17),
                       limits = c("HD", "CW"),
                       labels = c("Hot & Dry", "Cool & Wet")) + 
    theme(text = element_text(size = 20), 
          panel.grid.major = element_line(colour = "grey"))
} 
# End function -----------------------------------------------------------------

# plot_gene_counts_day(dds, gene = "MA_135205g0010")

gene_mod_list <- readRDS("gene_mod_list.rds")

sig_mod_colors <- c("darkgrey", "magenta", "darkred", "red")

for (color in sig_mod_colors) {

  plot_list <- vector(mode = "list", 
                      length = length(gene_mod_list[[color]]))
  
  for (gene in gene_mod_list[[color]]){
    g <- plot_gene_counts_day(dds, gene = gene[i])
    plot_list[[gene]] <- g
  }
  
  pdf(paste0(color, "_all_genes_counts_day.pdf"), 
      width = 10, height = 6)
  
  for(gene in gene_mod_list[[color]]){
    print(plot_list[[gene]])
  }
  dev.off()
}

# Tryna plot all on one --------------------------------------------------------

# Packages
require(tidyverse)
source(here::here("myscripts/functions.R"))

# List of files to load
setwd(here::here("myresults/gwas_env"))
files <- list.files()[str_which(list.files(), "res_")]
names <- str_replace(str_extract(files, "res.+"), ".rds", "")

# Read in .rds file and assign each object to the base name of the file
for (i in seq_along(files)) {
  assign(names[i], readRDS(files[i]))
}


# Bind together all the relevant res objects

obj <- objects()[c(str_which(objects(), "HD"), str_which(objects(), "CW"))]
dfs <- vector(mode = "list", length = length(obj))

for (i in seq_along(obj)) {
  dfs[[i]] <- as_tibble(get(obj[i]), 
                        rownames = "gene") %>%
    mutate(res = obj[i])
}


df_all <- bind_rows(dfs) %>%
  separate(res, 
           into = c("res", "climate", "treatment", "day", "ref0"), 
           sep = "_") %>%
  mutate(day = as.numeric(day)) %>%
  mutate(module = case_when(gene %in% gene_mod_list[["darkgrey"]] ~ "darkgrey",
                            gene %in% gene_mod_list[["darkred"]] ~ "darkred",
                            gene %in% gene_mod_list[["magenta"]] ~ "magenta",
                            gene %in% gene_mod_list[["red"]] ~ "red")) %>%
  filter(!is.na(module)) %>%
  add_row(log2FoldChange = rep(0, 3), 
          day = rep(0, 3), 
          treatment = c("C", "D", "H"),
          module = "darkgrey") %>%
  add_row(log2FoldChange = rep(0, 3), 
          day = rep(0, 3), 
          treatment = c("C", "D", "H"),
          module = "darkred") %>%
  add_row(log2FoldChange = rep(0, 3), 
          day = rep(0, 3), 
          treatment = c("C", "D", "H"),
          module = "magenta") %>%
  add_row(log2FoldChange = rep(0, 3), 
          day = rep(0, 3), 
          treatment = c("C", "D", "H"),
          module = "red")

df_avg <- df_all %>%
  dplyr::group_by(module, day, treatment) %>%
  dplyr::mutate(abs_lfc = abs(log2FoldChange)) %>%
  dplyr::summarize(avg_lfc = mean(abs_lfc),
                   sd_lfc = sd(abs_lfc))


# ------------------------------------------------------------------------------
# Function: plot_avg_lfc_time
# Description: Plot the average log-fold change over time for each treatment
# Inputs: data.frame and module character string  
# Outputs: ggplot

plot_avg_lfc_time <- function(data) {
  ggplot(data, 
         aes(x = day, y = avg_lfc, color = treatment)) +
    geom_point(size = 3, position = position_dodge(width = 0.5)) +
    geom_line(position = position_dodge(width = 0.5)) + 
    geom_errorbar(aes(ymin = avg_lfc - sd_lfc,
                      ymax = avg_lfc + sd_lfc),
                  width = .2,
                  position = position_dodge(width = 0.5)) +
    scale_color_manual(name = "Treatment",
                       values = c("grey50", "blue", "firebrick"),
                       limits = c("C", "H", "D"),
                       labels = c("Control", "Heat", "Heat & Drought")) + 
    geom_hline(yintercept = 0) +
    labs(x = "Day",
         y = "Absolute value log-fold change\nrelative to day 0") +
    ylim(c(-1, 4)) + 
    xlim(-.5,11) +
    theme_classic(base_size = 14) + 
    facet_wrap(~ module, ncol = 2)
} 
# End function -----------------------------------------------------------------

plot_avg_lfc_time(df_avg)
  


