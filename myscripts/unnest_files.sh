#! /bin/bash

# script to unnest the directories for the files for Emily's 110 transcriptomics samples into one folder with all files

# cd /Volumes/lockwood_lab/H202SC..../Rawdata/

mkdir /Volumes/lockwood_lab/emily_rnaseq

for file in /Volumes/lockwood_lab/H202SC19120714/Rawdata/*/*.gz

  do
    f=`basename ${file}`
    mv ${file} /Volumes/lockwood_lab/emily_rnaseq/${f}
  done
