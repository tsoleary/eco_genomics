# ------------------------------------------------------------------------------
# Over-representation analysis on P.rubens transcriptomics data
# April 22, 2020
# Katelynn Warner, Ben Camber, & TS O'Leary
# ------------------------------------------------------------------------------

# Load packages ----------------------------------------------------------------
require(dplyr)
require(ggplot2)
require(reshape2)
require(topGO)
require(plyr)

# Read in some table with gene ids, for example expressions
setwd(here::here("myresults/gwas_env"))

geneID2GO <- readMappings("Pabies1.0-gene_go_concat.txt")

geneNames <- names(geneID2GO)

# Load gene mod list
gene_mod_list <- readRDS("gene_mod_list.rds")

sig_mod_colors <- c("darkgrey", "magenta", "darkred", "red")

for (color in sig_mod_colors) {
  
  # Get the list of genes of interest
  geneList <- factor(as.integer(geneNames %in% gene_mod_list[[color]]))
  names(geneList) <- geneNames
  
  # Load in data with each gene list
  GOdata <- new("topGOdata", 
                ontology = "BP", 
                allGenes = geneList, 
                node = 10,
                annot = annFUN.gene2GO, 
                gene2GO = geneID2GO)
  
  # Run topGO with elimination test
  resultTopGO.elim <- runTest(GOdata, 
                              algorithm = "elim", 
                              statistic = "Fisher")
  
  allRes <- GenTable(GOdata, 
                     elimKS = resultTopGO.elim,
                     orderBy = "elimKS", 
                     topNodes = 1000) %>%
    filter(elimKS < 0.05 & Annotated <= 200)
  
  # Write and save topGO_results
  write.table(allRes, 
              file = paste(color, "topGO_results.txt", sep = "_"), 
              sep = "\t", 
              quote = FALSE, 
              col.names = TRUE, 
              row.names = FALSE)
}

# Load the ORA results from the WGCNA results ----------------------------------
setwd(here::here("myresults/gwas_env"))
list.files()[str_detect(list.files(), "topGO_results.txt")]

# List of files to load
setwd(here::here("myresults/gwas_env"))
files <- list.files()[str_detect(list.files(), "topGO_results.txt")]
names <- paste0("ora_", str_replace(files, "_topGO_results.txt", ""))

# Read in .rds file and assign each object to the base name of the file
for (i in seq_along(files)) {
  assign(names[i], read_delim(files[i], delim = "\t"))
}

