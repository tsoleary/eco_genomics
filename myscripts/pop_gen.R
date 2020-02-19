
# setwd

# load the _outFold.sfs file

SFS <- scan("BFA_outFold.sfs")

#
sumSFS <- sum(SFS)

# the first column of the site frequency spectrum is the monomorphic sites
pctPoly <- 100*(1 - (SFS[1]/sumSFS))


# read in the pestPG table
div <- read.table("BFA_folded_allstites.thetas.idx.pestPG")

colnames(div) <- c("window", "chrname", "wincenter", 
                   "tW", "tP", "tF", "tH", "tL", "tajD",
                   "fulif", "fuliD", "fayH", "zengsE", "numSites")

# tW / numSites
# tP / numSites
                
                