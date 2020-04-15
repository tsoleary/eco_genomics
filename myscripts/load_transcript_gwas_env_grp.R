# ------------------------------------------------------------------------------
# Load transcriptomics data related to the GWAS-Env group project
# April 15, 2020
# TS O'Leary
# ------------------------------------------------------------------------------

# Packages
require(tidyverse)

# List of files to load
setwd(here::here("myresults/gwas_env"))
res_files <- list.files()[str_which(list.files(), "res")]
res_names <- str_replace(str_extract(res_files, "res.+"), ".rds", "")

# Read in .rds file and assign each object to the base name of the file
for (i in seq_along(res_files)) {
  assign(res_names[i], readRDS(res_files[i]))
}
