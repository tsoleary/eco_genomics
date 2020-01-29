#! /bin/bash/

# make these results in myresults directory
cd ~/eco_genomics/myresults/

# make a new fastqc directory
mkdir fastqc


# run a for loop to do the fastqc on all the files within the populations
for file in /data/project_data/RS_ExomeSeq/fastq/edge_fastq/BFA*fastq.gz

do 

fastqc ${file} -o fastqc/

# close the for loop
done

