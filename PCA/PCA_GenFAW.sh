#!/bin/bash
#SBATCH -p cpu-dedicated
#SBATCH --account=dedicated-cpu@dgimi-eha
#SBATCH --mem=64G
#SBATCH --cpus-per-task=8

set -euo pipefail
# ==============================
# ENVIRONNEMENTS
# ==============================
source /home/durandk/miniconda3/etc/profile.d/conda.sh
conda activate plink2
# ==============================
# PCA #--geno 0.05: eliminates SNPs with more than 5% missing, retaining those with at least 95% of data present.

# ==============================
plink2 --vcf  /storage/simple/projects/faw_adaptation/Data_Backup/Merged_vcf/2025_GenFAW600/2025_GenFAW600/VCF_pruned_for_PCA_Admixture/WholeGenome_biallelic_max5miss_pruned.vcf.gz  \
	  --chr-set 29 \
	  --allow-extra-chr \
	  --pca \
	  --out WholeGenome_biallelic_max5miss_pruned_PCA



#############################################
# PCA autosomes
# #############################################

plink2 --vcf  /storage/simple/projects/faw_adaptation/Data_Backup/Merged_vcf/2025_GenFAW600/2025_GenFAW600/VCF_pruned_for_PCA_Admixture/Autosome_biallelic_max5miss_pruned.vcf.gz \
	  --chr-set 29 \
	  --allow-extra-chr \
	  --pca \
	  --out Autosome_biallelic_max5miss_pruned_PCA
	 
#############################################
# PCA Z chromosome
# #############################################

plink2 --vcf  /storage/simple/projects/faw_adaptation/Data_Backup/Merged_vcf/2025_GenFAW600/2025_GenFAW600/VCF_pruned_for_PCA_Admixture/Z_biallelic_max5miss_pruned.vcf.gz \
	  --allow-extra-chr \
	  --pca \
	  --out Z_biallelic_max5miss_pruned_PCA
