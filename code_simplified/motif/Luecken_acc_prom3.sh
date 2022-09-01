#!/bin/bash

#$ -cwd
#$ -pe smp 1
#$ -l m_mem_free=64G
#$ -l h_rt=25:00:00
#$ -N motif

source ~/.bash_profile
path_to_file=/fast/AG_Ohler/hyliou/mtproject/processed/Luecken/acc_prom1
files=$path_to_file/*
out=/fast/AG_Ohler/hyliou/mtproject/Luecken/acc_prom_motif

mkdir -p $out

for file in $files
do
    bsfile=$(basename $file)
    mkdir -p ${out}/${bsfile%.bed}
    findMotifsGenome.pl ${path_to_file}/${bsfile} hg38 ${out}/${bsfile%.bed} -bg /fast/AG_Ohler/hyliou/mtproject/processed/Luecken_ATAC_prom.bed
done