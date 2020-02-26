#! /bin/bash
# this is a script to rename a bunch of files I accidentially messed up the naming for

# path to the folder that they are in
output="/data/project_data/RS_ExomeSeq/mapping"

# rename the .sorted.rmdup.dup.bam files to .sorted.rmdup.bam

for file in ${output}/BWA/*.sorted.rmdup.dup.bam

do
  f=${file/.sorted.rmdup.dup.bam/.sorted.rmdup.bam}
  mv ${file} ${f}
done

# rename the .sorted.rmdup.dup.bam.bai files to .sorted.rmdup.bam.bai
for file in ${output}/BWA/*.sorted.rmdup.dup.bam.bai

do
  f=${file/.sorted.rmdup.dup.bam.bai/.sorted.rmdup.bam.bai}
  mv ${file} ${f}
done
