#!/bin/bash

bismark --bowtie2 --multicore 1 \
    --genome /data/project_data/epigenetics/reference_genome \
    --output_dir /data/project_data/epigenetics/aligned_output \
    -1 /data/project_data/epigenetics/trimmed_fastq/AH_F25_1_1.fq.gz \
    -2 /data/project_data/epigenetics/trimmed_fastq/AH_F25_1_2.fq.gz \
    --rg_tag --rg_id AH_F25_1 --rg_sample AH_F25_1 --gzip --local --maxins 1000
