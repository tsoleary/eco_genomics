#! /bin/bash

# this script is to run the read mapping using the bwa program

# -t 1 is that we are using only one thread
# -M is a special flag if its mapping is split across more than one contig
# this just allows it to map across contigs
# ${ref} variable ref that we will define as the reference genome
# then it will redirect the output file .sam to the output dir we made
# > redirects the output to save it to a file, otherwise bash would just vomit
# the words on the screen

ref="/data/project_data/RS_ExomeSeq/ReferenceGenomes/Pabies1.0-genome_reduced.fa"

# a loop to map each individual within my population
# ${input} is the path to the

for forward in ${input}*_R1.cl.pd.fq

do
  # define the reverse read based on that forward read
  # then define the forward read as the whole text
  # and replace it with nothing so that you only have the base part of the name
  # having the slash inside the {} means basically find and replace... so the
  # first part of the slash is interpreted as the find and the second part is
  # what to replace it with
  reverse=${forward/_R1.cl.pd.fq/_R2.cl.pd.fq}
  f=${forward/_R1.cl.pd.fq/}
  # basename is a function in base that creates the basename that will be used
  # for naming the rest of the files
  name=`basename ${f}`
  bwa mem -t 1 -M ${ref} ${forward} ${reverse} > ${output}/BWA/${name}.sam
done
