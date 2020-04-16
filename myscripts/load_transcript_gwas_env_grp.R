# ------------------------------------------------------------------------------
# Load transcriptomics data related to the GWAS-Env group project
# April 15, 2020
# TS O'Leary
# ------------------------------------------------------------------------------

# Packages
require(tidyverse)
source(here::here("myscripts/functions.R"))

# List of files to load
setwd(here::here("myresults/gwas_env"))
files <- list.files()[str_which(list.files(), "res")]
names <- str_replace(str_extract(files, "res.+"), ".rds", "")

# Read in .rds file and assign each object to the base name of the file
for (i in seq_along(files)) {
  assign(names[i], readRDS(files[i]))
}


# DEG lists --------------------------------------------------------------------

fdr <- 0.05

names <- objects()[str_which(objects(), "res")]

for (i in seq_along(names)) {
  assign(str_replace(names[i], "res", "deg"), 
         extract_deg(get(names[i]), fdr = fdr))
}


intersect(deg_5_0[["all"]], deg_10_0[["all"]])

intersect(deg_5_0[["up"]], deg_10_0[["up"]])

intersect(deg_5_0[["down"]], deg_10_0[["down"]])






