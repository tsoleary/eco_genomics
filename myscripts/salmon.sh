#! /bin/bash/

# run salmon to do the mapping and the counting of the cleanreads

# salmon quant -i transcripts_index -l <LIBTYPE> -r reads.fq --validateMappings -o transcripts_quant

cd /data/project_data/RS_RNASeq/fastq/cleanreads/

for file in CAM_02_H_*.cl.fq

do

  salmon quant \
  -i /data/project_data/RS_RNASeq/ReferenceTranscriptome/Pabies_cds_index \
  -l A \
  -r ${file} \
  --validateMappings \
  -o /data/project_data/RS_RNASeq/salmon/cleanedreads/${file}

done
