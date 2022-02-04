#!/bin/bash

gene1=$1          # 1st positional argument in BASH
gene2=$2
output_name=$3

GENE_FASTA_DIR='/home/andrew/crisporWebsite/sampleFiles'

gene1='musculusgenes/Mus_musculus_ENSMUST00000108657_4_sequence.fa'     # Location of the p53 FASTA sequence for Musculus
gene2='spretusgenes/Mus_spretus_MGP_SPRETEiJ_T0028933_1_sequence.fa'    #  " " for Spretus
output_name='p53'

species1='ens90MusMus'      # Mus Musculus
species2='GCA_001624865.1'  # Mus Spretus
output_suffix='.txt'        # suffix to append to output. Could also be .tsv

touch $GENE_FASTA_DIR/$output_name+$species1+$version+$output_suffix
touch $GENE_FASTA_DIR/$output_name+$species2+$version+$output_suffix
touch $output_name+$species1+$species2+$version

versions=('-m' '-control')  # Finds unique guides, then homozygous KO guides
for version in ${versions[@]}; do 
python3.6 findingsgRNAs.py $version \
  $GENE_FASTA_DIR/$gene1 \
  $species1 \
  $GENE_FASTA_DIR/$output_name+$species1+$version+$output_suffix \
  $GENE_FASTA_DIR/$gene2 \
  $species2 \
  $GENE_FASTA_DIR/$output_name+$species2+$version+$output_suffix \
  $output_name+$species1+$species2+$version
done
