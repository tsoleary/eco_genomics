# ------------------------------------------------------------------------------
# Load epigenetics data for epigenetics write up
# April 08, 2020
# TS O'Leary
# ------------------------------------------------------------------------------

require(methylKit)
require(tidyverse)

# Set directory
dir <- here::here("myresults/epigenetics")


# Mapping rate info ------------------------------------------------------------

# map_rate <- read_csv("alignment_rate.csv") %>%
#   separate(Sample, into = c("Group", "Gen", "Rep", "One"), sep = "_",
#            remove = FALSE) %>% filter(Gen != "F00")
# 
# ggplot(map_rate, aes(x = Group, y = Mapping_percentage, fill = Group)) +
#   geom_jitter(width = 0.1, pch = 21, size = 3, color = "grey50", alpha = 0.8) +
#   ylim(c(25, 65)) +
#   ylab("Mapping rate (%)") +
#   xlab("") +
#   scale_fill_manual(name = "Treatment",
#                     values = wesanderson::wes_palette("GrandBudapest1"),
#                     labels = c("Control", 
#                                expression("CO"[2]), 
#                                "Heat", 
#                                expression("Heat & CO"[2]))) +
#   scale_x_discrete(limits = c("AA", "AH", "HA", "HH"),
#                    labels = c("Control", 
#                               expression("CO"[2]), 
#                               "Heat", 
#                               expression("Heat & CO"[2]))) +
#   theme_classic(base_size = 14) +
#   theme(legend.text.align = 0)
# 
# map_rate %>%
#   summarize(avg_mapping_rate = mean(Mapping_percentage))
# 
# map_rate %>%
#   group_by(Group, Gen) %>%
#   summarize(avg_mapping_rate = mean(Mapping_percentage))

# Load methyl ------------------------------------------------------------------
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

# Subset meth with comparisons -------------------------------------------------
meth_sub_HA <- reorganize(meth,
                          sample.ids = c("AA_F25_1", 
                                         "AA_F25_2", 
                                         "AA_F25_3", 
                                         "AA_F25_4",
                                         "HA_F25_1",
                                         "HA_F25_2",
                                         "HA_F25_3",
                                         "HA_F25_4"),
                          treatment = c(0,0,0,0,1,1,1,1),
                          save.db = FALSE)

meth_sub_HH <- reorganize(meth,
                          sample.ids = c("AA_F25_1", 
                                         "AA_F25_2", 
                                         "AA_F25_3", 
                                         "AA_F25_4",
                                         "HH_F25_1",
                                         "HH_F25_2",
                                         "HH_F25_3",
                                         "HH_F25_4"),
                          treatment = c(0,0,0,0,1,1,1,1),
                          save.db = FALSE)


# Calculate differential methylation -------------------------------------------
myDiff_HA <- calculateDiffMeth(meth_sub_HA,
                               overdispersion = "MN",
                               mc.cores = 1,
                               suffix = "AA_HH",
                               adjust = "qvalue",
                               test = "Chisq")

myDiff_HH <- calculateDiffMeth(meth_sub_HH,
                               overdispersion = "MN",
                               mc.cores = 1,
                               suffix = "AA_HH",
                               adjust = "qvalue",
                               test = "Chisq")

# Get all differentially methylated bases
myDiff_HA <- getMethylDiff(myDiff_HA, qvalue = 0.05, difference = 10)
myDiff_HH <- getMethylDiff(myDiff_HH, qvalue = 0.05, difference = 10)

# Percent methylation 
pm_HA <- percMethylation(meth_sub_HA)
# make a dataframe with snp id's, methylation, etc.
sig.in_HA <- as.numeric(row.names(myDiff_HA))
pm.sig_HA <- pm_HA[sig.in_HA, ]

# add snp, chr, start, stop
din_HA <- getData(myDiff_HA)[,1:3]
df.out_HA <- cbind(paste(getData(myDiff_HA)$chr, 
                         getData(myDiff_HA)$start, 
                         sep =  ":"), 
                   din_HA, pm.sig_HA)
colnames(df.out_HA) <- c("snp", colnames(din_HA), 
                         colnames(df.out_HA[5:ncol(df.out_HA)]))
df.out_HA <- (cbind(df.out_HA, getData(myDiff_HA)[, 5:7]))


# # write bed file for intersection with genome annotation
# setwd(here::here("myresults/epigenetics"))
# write.table(file = "HA_diffmeth.bed",
#             data.frame(chr= df.out_HA$chr, 
#                        start = df.out_HA$start, 
#                        end = df.out_HA$end),
#             row.names = FALSE, 
#             col.names = FALSE, 
#             quote = FALSE, 
#             sep = "\t")


# Percent methylation 
pm_HH <- percMethylation(meth_sub_HH)
# make a dataframe with snp id's, methylation, etc.
sig.in_HH <- as.numeric(row.names(myDiff_HH))
pm.sig_HH <- pm_HH[sig.in_HH, ]

# add snp, chr, start, stop
din_HH <- getData(myDiff_HH)[,1:3]
df.out_HH <- cbind(paste(getData(myDiff_HH)$chr, 
                         getData(myDiff_HH)$start, 
                         sep =  ":"), 
                   din_HH, pm.sig_HH)
colnames(df.out_HH) <- c("snp", colnames(din_HH), 
                         colnames(df.out_HH[5:ncol(df.out_HH)]))
df.out_HH <- (cbind(df.out_HH, getData(myDiff_HH)[, 5:7]))


# # write bed file for intersection with genome annotation
# setwd(here::here("myresults/epigenetics"))
# write.table(file = "HH_diffmeth.bed",
#             data.frame(chr= df.out_HH$chr, 
#                        start = df.out_HH$start, 
#                        end = df.out_HH$end),
#             row.names = FALSE, 
#             col.names = FALSE, 
#             quote = FALSE, 
#             sep = "\t")
# 
# # Run epigenetics_link2gene.sh on the server to get the hits.bed files

# Read the tables back into R --------------------------------------------------
setwd(here::here("myresults/epigenetics"))
hits_HA <- read_delim("HA_hits.bed", 
                      delim = "\t", 
                      col_names = c("SNP", "Start", "Stop", "ScaffoldName", 
                                    "FromPosition", "ToPosition",  
                                    "Sense", "TranscriptName", "TranscriptPath",  
                                    "GeneAccession", "GeneName", "GeneAltNames",  
                                    "GenePfam", "GeneGo", "CellularComponent", 
                                    "MolecularFunction", "BiologicalProcess",
                                    "PositionFromGene"))

hits_HH <- read_delim("HH_hits.bed", 
                      delim = "\t", 
                      col_names = c("SNP", "Start", "Stop", "ScaffoldName", 
                                    "FromPosition", "ToPosition",  
                                    "Sense", "TranscriptName", "TranscriptPath",  
                                    "GeneAccession", "GeneName", "GeneAltNames",  
                                    "GenePfam", "GeneGo", "CellularComponent", 
                                    "MolecularFunction", "BiologicalProcess",
                                    "PositionFromGene"))

# Create dfs for plots ---------------------------------------------------------

# # Individual gene plot -----
# HA_ind_plot <- df.out_HA %>%
#   select(snp, contains("_")) %>%
#   pivot_longer(-snp, names_to = "sample", values_to = "methylation") %>%
#   separate(sample, into = c("Group", "Gen", "Rep"), sep = "_")
# 
# HA_ind_plot %>% 
#   filter(snp == "LS043685.1:10221") %>%
#   ggplot(aes(x = Group, y = methylation, fill = Group)) +
#   stat_summary(color = "grey50", fun.data = "mean_se") +
#   geom_jitter(width = 0.1, size = 3, pch = 21, alpha = 0.8) + 
#   ylab("Methylation (%)") +
#   xlab("Treatment") +
#   scale_fill_manual(values = c("grey50", "pink")) +
#   scale_x_discrete(limits = c("AA", "HA"),
#                    labels = c("Control", 
#                               "Heat")) +
#   theme_classic(base_size = 18) +
#   theme(legend.position = "none")


# Percent methylation plot df -----
pm <- percMethylation(meth) %>%
  as_tibble() %>%
  pivot_longer(everything(), values_to = "methylation",
               names_to = "sample") %>%
  separate(sample, into = c("Group", "Gen", "Rep"), sep = "_") %>%
  filter(Gen != "F00" & Group != "AH") %>%
  group_by(Group, Rep) %>%
  summarize(avg_meth = mean (methylation))

# Diff meth plot histogram -----

# Calc numbers hypo and hypermethylated
# df.out_HA %>%
#   group_by(meth.diff < 0) %>%
#   count()
# 
# df.out_HH %>%
#   group_by(meth.diff < 0) %>%
#   count()

df.out_HA$Treatment <- "HA"
df.out_HH$Treatment <- "HH"

hist.plot <- bind_rows(df.out_HA, df.out_HH)

# KS Test on the distributions
KS_pval <- ks.test(hist.plot$meth.diff[hist.plot$Treatment == "HA"],
                   hist.plot$meth.diff[hist.plot$Treatment == "HH"])$p.value

# Differentially methylated gene lists -----------------------------------------
HA_dmg <- hits_HA %>%
  filter(GeneName != ".") %>%
  distinct(GeneName, .keep_all = TRUE)

HH_dmg <- hits_HH %>%
  filter(GeneName != ".") %>%
  distinct(GeneName, .keep_all = TRUE)
