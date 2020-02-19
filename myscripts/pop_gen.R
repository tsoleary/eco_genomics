
require(tidyverse)

# setwd
setwd(here::here("myresults/ANGSD"))

# load the _outFold.sfs file
SFS <- scan("BFA_outFold.sfs")

# this is the 
sumSFS <- sum(SFS)

# the first column of the site frequency spectrum is the monomorphic sites
pctPoly <- 100*(1 - (SFS[1]/sumSFS))


# read in the pestPG table -------
div <- read.table("BFA_folded_allsites.thetas.idx.pestPG")

colnames(div) <- c("window", "chrname", "wincenter", 
                   "tW", "tP", "tF", "tH", "tL", "tajD",
                   "fulif", "fuliD", "fayH", "zengsE", "numSites")
# tW is thetaW or waterson's theta
# tP is thetaPi or theta pi



# tW / numSites # tP / numSites
div <- div %>%
  mutate(tWperSite = tW / numSites) %>%
  mutate(tPperSite = tP / numSites)

# take the mean for the tWperSite, tPperSite, and tajD
# or just get it from the mean in the specific column using summary(div) 
div %>%
  summarize(avgtWperSite = mean(tWperSite),
            avgtPperSite = mean(tPperSite),
            avgtajD = mean(tajD))


# visualize the distribution of values across the sites
ggplot(SFS) +
  geom_bar(aes(x = SFS))

barplot(SFS[-1]) # remove the monomorphic site to visualize

# plot the tWperSite distribution
p1 <- ggplot(div) +
  geom_histogram(aes(x = tWperSite), color = "grey", alpha = 0.75) +
  theme_classic()

# plot the tPperSite distribution
p2 <- ggplot(div) +
  geom_histogram(aes(x = tPperSite), color = "grey", alpha = 0.75) +
  theme_classic()

# combine these plots together
p_1_2 <- cowplot::plot_grid(p1, p2, nrow = 1, labels = "AUTO")

# plot the tajD distribution
p3 <- ggplot(div) +
  geom_histogram(aes(x = tajD), color = "grey", alpha = 0.75) +
  theme_classic()

# combine all plots in one view
cowplot::plot_grid(p_1_2, p3, nrow = 2, rel_widths = c(1,2), labels = c("", "C"))

  


                
                