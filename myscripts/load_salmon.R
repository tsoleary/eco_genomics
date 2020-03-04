# script to load the quant.sf files into R and save as a single gene count table

# These commands need to be run in R on the server in command-line
# open R with the command:    R



library(tximportData)
library(tximport)

#locate the directory containing the files. 
dir <- "/data/project_data/RS_RNASeq/salmon/allmapping/"
list.files(dir)

# read in table with sample ids
samples <- read.table("/data/project_data/RS_RNASeq/salmon/RS_samples.txt", header=TRUE)

# now point to quant files
all_files <- file.path(dir, samples$sample, "quant.sf")
names(all_files) <- samples$sample

# what would be used if linked transcripts to genes
#txi <- tximport(files, type = "salmon", tx2gene = tx2gene)
# to be able to run without tx2gene
txi <- tximport(all_files, type = "salmon", txOut=TRUE)  
names(txi)

head(txi$counts)

countsMatrix <- txi$counts
dim(countsMatrix)
#[1] 66069    76

# To write out
write.table(countsMatrix, file = "RS_countsMatrix.txt", col.names = T, row.names = T, quote = F) 