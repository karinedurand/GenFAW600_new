#!/bin/bash
#SBATCH -p cpu-dedicated
#SBATCH --account=dedicated-cpu@dgimi-eha
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

#############################################
# 1. Get VCF file from vcf.gz
#############################################

# gunzip  WholeGenome_biallelic_max5missing_pruned.vcf.gz
#############################################
# 2. Convert PED to GENO for sNMF
#############################################

source /home/durandk/miniconda3/etc/profile.d/conda.sh
conda activate R

Rscript ped2geno.R 


