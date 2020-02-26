#! /bin/bash/

# make these results in myresults directory
cd ~/eco_genomics/myresults/fastqc/

# make a new fastqc directory
mkdir RNASeq_trim

# run a for loop to do the fastqc on all the files for my CAM pop
for file in /data/project_data/RS_RNASeq/fastq/CAM_02_H_*.cl.fq

do

  fastqc ${file} -o RNASeq/

# close the for loop
done

# my ESC pops ------------------------------------------------------------------
for file in /data/project_data/RS_RNASeq/fastq/ESC_01_C_*.cl.fq

do

  fastqc ${file} -o RNASeq/

# close the for loop
done
