#! /bin/bash/

# script to visualize the quality of the trimmed reads

# make these results in myresults directory
cd ~/eco_genomics/myresults/

# make a new fastqc directory
mkdir trim_fastqc

# run a for loop to do the fastqc on all the files within the populations
for file in /data/project_data/RS_ExomeSeq/fastq/edge_fastq/pairedcleanreads/BFA*cl.pd.fq

do 

fastqc ${file} -o trim_fastqc/

# close the for loop
done

