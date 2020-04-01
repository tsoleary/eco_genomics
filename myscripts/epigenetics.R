# ------------------------------------------------------------------------------
# Epigenetics script for EcoGenomics Tutorial
# April 01, 2020
# TS O'Leary
# ------------------------------------------------------------------------------

library(methylKit)
library(tidyverse)
library(ggplot2)
library(pheatmap)

# Read in the raw methylation calls with methylkit -----------------------------

# Set directory
dir <- here::here("myresults/epigenetics")

# Read in the sample IDs
samples <- read.table(paste0(dir, "/", "sample_id.txt"), 
                      header = FALSE)
# Point to coverage files
files <- file.path(dir, samples$V1)
all(file.exists(files))

# Convert to list
file.list <- as.list(files)
# Get the names only for naming our samples
nmlist <- as.list(gsub("_1_bismark_bt2_pe.bismark.cov.gz", "", samples$V1))

# use methRead to read in the coverage files
myobj <- methRead(location = file.list,
                  sample.id = nmlist,
                  assembly = "atonsa",
                  dbtype = "tabix",
                  context = "CpG",
                  resolution = "base",
                  mincov = 20,
                  treatment = c(0, 0, 0, 0,
                                1, 1, 1, 1,
                                2, 2, 2, 2,
                                3, 3, 3, 3,
                                4, 4, 4, 4),
                  pipeline = "bismarkCoverage",
                  dbdir = dir)

# Visualize coverage and filter ------------------------------------------------

# We can look at the coverage for individual samples with getCoverageStats()
getCoverageStats(myobj[[1]], plot = TRUE)

# Plot all of our samples at once to compare.


# Filter samples by depth
filtered.myobj <- filterByCoverage(myobj, 
                                   lo.count = 20,
                                   lo.perc = NULL,
                                   hi.count = NULL,
                                   hi.perc = 97.5,
                                   dbdir = dir)

# Merge samples ----------------------------------------------------------------
# Merge all the samples. We will require sites to be present in each sample
# Note this takes a couple of hours to run
meth <- unite(filtered.myobj, mc.core = 3, suffix = united, dbdir = dir)


# Load a previously created database -------------------------------------------
# So Reid created this and we are going to move on from here
meth <- methylKit:::readMethylBaseDB(
  dbpath = paste0(dir, "/methylBase_united.txt.bgz"),
  dbtype = "tabix",
  sample.id =   unlist(nmlist),
  assembly = "atonsa", # this is just a string. no actual database
  context = "CpG",
  resolution = "base",
  treatment = c(0,0,0,0,
                1,1,1,1,
                2,2,2,2,
                3,3,3,3,
                4,4,4,4),
  destrand = FALSE)

# Methylation statistics across samples ----------------------------------------

# percMethylation() calculates the percent methylation for each site and sample
pm <- percMethylation(meth)

#plot methylation histograms
ggplot(gather(as.data.frame(pm)), aes(value)) + 
  geom_histogram(bins = 10, color = "black", fill = "grey") + 
  facet_wrap(~ key, ncol = 4)

# calculate and plot mean methylation
sp.means <- colMeans(pm)

p.df <- data.frame(sample = names(sp.means),
                   group = substr(names(sp.means), 1, 6),
                   methylation = sp.means)

ggplot(p.df, aes(x = group, y = methylation, color = group)) + 
  stat_summary(color = "black") + 
  geom_jitter(width = 0.1, size = 3) +
  theme_classic()

# Summarize variation ----------------------------------------------------------

# Sample clustering
clusterSamples(meth, dist = "correlation", method = "ward.D", plot = TRUE)

# PCA
PCASamples()


# Find differentially methylated sites between two groups ----------------------

# Subset groups
meth_sub <- reorganize(meth,
                       sample.ids = c("AA_F00_1","AA_F00_2","AA_F00_3", "AA_F00_4",
                                      "HH_F25_1","HH_F25_2","HH_F25_3","HH_F25_4"),
                       treatment = c(0,0,0,0,1,1,1,1),
                       save.db = FALSE)

# Calculate differential methylation
myDiff <- calculateDiffMeth(meth_sub,
                            overdispersion = "MN",
                            mc.cores = 1,
                            suffix = "AA_HH",
                            adjust = "qvalue",
                            test = "Chisq")

# Get all differentially methylated bases
myDiff <- getMethylDiff(myDiff, qvalue = 0.05, difference = 10)

# We can visualize the changes in methylation frequencies quickly.
hist(getData(myDiff)$meth.diff)

# Hyper-methylated bases

# Hypo-methylated bases


# Plots of differentially methylated groups ------------------------------------
# Heatmaps -----
pm <- percMethylation(meth_sub)
# make a dataframe with snp id's, methylation, etc.
sig.in <- as.numeric(row.names(myDiff))
pm.sig <- pm[sig.in, ]

# add snp, chr, start, stop
din <- getData(myDiff)[,1:3]
df.out <- cbind(paste(getData(myDiff)$chr, 
                      getData(myDiff)$start, 
                      sep= ":"), 
                din, pm.sig)
colnames(df.out) <- c("snp", colnames(din), colnames(df.out[5:ncol(df.out)]))
df.out <- (cbind(df.out,getData(myDiff)[, 5:7]))

pheatmap(pm.sig, show_rownames = FALSE)

# We can also normalize 



# Methylation of specific snps -------------------------------------------------

# convert data frame to long form

# write bed file for intersection with genome annotation

