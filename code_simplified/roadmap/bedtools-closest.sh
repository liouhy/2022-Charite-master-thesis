#!/bin/bash
# find closest 

MY_PATH=/fast/AG_Ohler/hyliou/mtproject/

source ~/miniconda3/etc/profile.d/conda.sh
conda activate mtproject

for file in ${MY_PATH}raw/roadmap/*.gappedPeak.gz; do

    basename=$(basename $file)

    # with 10X PBMC ATAC data
    bedtools closest -a <(bedtools sort -i $file) -b ${MY_PATH}processed/10X_multiome/ATAC.bed -d -t first > ${MY_PATH}processed/roadmap/${basename%.gappedPeak.gz}.bed

done



