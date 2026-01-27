#!/bin/bash
#SBATCH -p cpu-dedicated
#SBATCH --account=dedicated-cpu@dgimi-eha
#SBATCH --cpus-per-task=8
#SBATCH --mem=64G

#############################################
# 1. Get VCF file from vcf.gz
#############################################


cp /storage/simple/projects/faw_adaptation/Data_Backup/Merged_vcf/2025_GenFAW600/2025_GenFAW600/VCF_pruned_for_PCA_Admixture/GenFAW_max5miss_pruned.vcf.gz  .
gunzip GenFAW_max5miss_pruned.vcf.gz  


#############################################
# 2. Convert PED to GENO for sNMF
#############################################

source /home/durandk/miniconda3/etc/profile.d/conda.sh
conda activate R

Rscript ped2geno.R 


