#!/bin/bash   

# this is a script for trimomatic for paired end reads of the data

# go to the directory that has the fastqc data
cd /data/project_data/RS_ExomeSeq/fastq/edge_fastq  

# we are going to start with the R1 forward read files 
for R1 in BFA*R1_fastq.gz  

# between the do and the done is where the action happens
do 
 
# there are some files that have pairs and some reads that don't have pairs so we want to have trimmomatic do it at the same time

        # define a second variable that is equal to the R1, 
        # find this bit of the name in R1 and replace it with the second with with the 
	R2=${R1/_R1_fastq.gz/_R2_fastq.gz}
	f=${R1/_R1_fastq.gz/}
	name=`basename ${f}`

# this is a java program and that program is located in that popgen directory
# the slashes allow you to write it on multiple lines but it will 
# interpret it as allon one line

# then there are a bunch of different directories that are defined to store the paired end clean or the unpaired for each of the R1 and R2

# trim of leading and trailing sequence that are less than 20
# set a window size of 6 and if the average is less than 20 ditch it
# minimum length is set to 35 bp

# then close the loop with the done


	java -classpath /data/popgen/Trimmomatic-0.33/trimmomatic-0.33.jar org.usadellab.trimmomatic.TrimmomaticPE \
        -threads 1 \
        -phred33 \
         "$R1" \
         "$R2" \
         /data/project_data/RS_ExomeSeq/fastq/edge_fastq/pairedcleanreads/${name}_R1.cl.pd.fq \
         /data/project_data/RS_ExomeSeq/fastq/edge_fastq/unpairedcleanreads/${name}_R1.cl.un.fq \
         /data/project_data/RS_ExomeSeq/fastq/edge_fastq/pairedcleanreads/${name}_R2.cl.pd.fq \
         /data/project_data/RS_ExomeSeq/fastq/edge_fastq/unpairedcleanreads/${name}_R2.cl.un.fq \
        ILLUMINACLIP:/data/popgen/Trimmomatic-0.33/adapters/TruSeq3-PE.fa:2:30:10 \
        LEADING:20 \
        TRAILING:20 \
        SLIDINGWINDOW:6:20 \
        MINLEN:35 
 
done 

