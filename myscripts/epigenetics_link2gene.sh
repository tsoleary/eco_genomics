#!/bin/bash

/data/popgen/bedtools2/bin/bedtools closest -a ~/eco_genomics/diffmeth.bed \
      -b /data/project_data/epigenetics/GCA_900241095.1_Aton1.0_genomic.fa_annotation_table.bed \
      -D b | \
      awk '!($10=="-")' > hits.bed

# the annotation file here has the format: ScaffoldName  FromPosition  ToPosition  Sense TranscriptName  TranscriptPath  GeneAccession GeneName  GeneAltNames  GenePfam  GeneGo  CellularComponent MolecularFunction BiologicalProcess

# note that the hits.bed file will paste the diffmeth.bed file before the annotation table. So the first 3 columns are fom diffmeth.bed, then next 8 from the annotation table.
