# ------------------------------------------------------------------------------
# EcoGenomics Functions
# TS O'Leary
# ------------------------------------------------------------------------------


# ------------------------------------------------------------------------------
# Function: extract_deg
# Description: Extract differentially expressed genes vector at fdr cutoff
# Inputs: results S4 object class of DESeq2, fdr cut off
# Outputs: character vector of degs

extract_deg <- function(res, fdr = 0.05) {
  deg_all <- rownames(subset(res, res$padj < fdr))
  deg_up <- rownames(subset(res, res$padj < fdr & res$log2FoldChange > 0))
  deg_down <- rownames(subset(res, res$padj < fdr & res$log2FoldChange < 0))
  
  list("all" = deg_all, "up" = deg_up, "down" = deg_down)
} 
# End function -----------------------------------------------------------------

