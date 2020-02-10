#! /bin/bash
# this is where out output sam files are going to get converted into binary
# format (bam)

# then we are going to sort the bam files, remove PCR duplicates, and index them

# first convert sam to bam and then sort

for f in ${output}/BWA/${mypop}*.sam

do

  out=${f/.sam/}
  # this converts it to the bam
  sambamba-0.7.1-linux-static view -S --format=bam ${f} -o ${out}.bam
  # now we need to sort the bam files so it can quickly find the function
  samtools sort ${out}.bam -o ${out}.sorted.bam

done

# next remove the pcr duplicates from our bam files

for file in ${output}/BWA/${mypop}*.sorted.bam

do

  f=${file/.sorted.bam/}
  # this is going to mark duplicates if they have the exact same coordinates, because they are randomly fragmented so if they star
  # -r removes all dups except for one
  # no dash o -o needed for the output
  # just calling the command on the server without any of the options gives you the help information
  sambamba-0.7.1-linux-static markdup -r -t 1 ${file} ${f}.sorted.rmdup.bam

done

# finish by indexing the files

for file in ${output}/BWA/${mypop}*sorted.rmdup.bam

do

  samtools index ${file}

done
