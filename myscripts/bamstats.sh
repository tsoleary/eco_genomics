#! /bin/bash

# set repo

myrepo="/users/t/s/tsoleary/eco_genomics"

mypop="BFA"

output="/data/project_data/RS_ExomeSeq/mapping"

# this echo is creating a header for the file and saving it on the flagstats.txt files
echo "Num.reads R1 R2 Paired MateMapped Singletons MateMappedDiffChr" > ${myrepo}/myresults/${mypop}.flagstats.txt

for file in ${output}/BWA/${mypop}*sorted.rmdup.bam
  do
    f=${file/.sorted.rmdup.bam/}
    name=`basename ${f}`
    # the double >> makes it so that it appends it and doesn't just rewrite the file every time
    echo ${name} >> ${myrepo}/myresults/${mypop}.names.txt
    # NR is the number of rows, row 6 through 12 have the meaningful results in them
    samtools flagstat ${file} | awk 'NR>=6&&NR<=12 {print $1}' | column -x
  done >> ${myrepo}/myresults/${mypop}.flagstats.txt

# calculate the depth of coverage from our bam files
for file in ${output}/BWA/${mypop}*sorted.rmdup.bam
  do
    samtools depth ${file} | awk '{sum+=$3} END {print sum/NR}'
  done >> ${myrepo}/myresults/${mypop}.coverage.txt
