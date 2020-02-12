#! /bin/bash

myrepo="/users/t/s/tsoleary/eco_genomics"

mypop="BFA"

mkdir ${myrepo}/myresults/ANGSD

output="${myrepo}/myresults/ANGSD"

# make a list of all of the files in this repo that match this and then save it as its own list file in the output file
ls /data/project_data/RS_ExomeSeq/mapping/BWA/${mypop}_*sorted.rm*.bam >${output}/${mypop}_bam.list

# we need to define where the reference genome is
REF="/data/project_data/RS_ExomeSeq/ReferenceGenomes/Pabies1.0-genome_reduced.fa"

# the back slashes tell bash to interpret it all on one line, just helps visuallize it. make sure there is nothing trailing the slash it needs to be a hard return.

# Estimating GL's and allele frequencies for all sites with ANGSD

ANGSD -b ${output}/${mypop}_bam.list \
-ref ${REF} -anc ${REF} \
-out ${output}/${mypop}_allsites \
-nThreads 1 \
-remove_bads 1 \
-C 50 \
-baq 1 \
-minMapQ 20 \
-minQ 20 \
-setMinDepth 3 \
-minInd 2 \
-setMinDepthInd 1 \
-setMaxDepthInd 17 \
-skipTriallelic 1 \
-GL 1 \
-doCounts 1 \
-doMajorMinor 1 \
-doMaf 1 \
-doSaf 1 \
-doHWE 1 \
# -SNP_pval 1e-6
