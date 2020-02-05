#! /bin/bash

# we will use this as a wrappeer to run our different scripts

# path to my repo on the server
myrepo="/users/t/s/tsoleary/eco_genomics"

# define the population
mypop="BFA"

# directory to our cleaned and paired paired
input="/data/project_data/RS_ExomeSeq/fastq/edge_fastq/pairedcleanreads/${mypop}"

# define an output directory in the common space not your own repo
output="/data/project_data/RS_ExomeSeq/mapping"


# call the mapping.sh and process_bam.sh files, you could also add the stats but for the sake of time we will do that later

source ./mapping.sh

source ./process_bam.sh
