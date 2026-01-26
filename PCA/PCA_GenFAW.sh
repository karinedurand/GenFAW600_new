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
# plink2 --vcf  /storage/simple/users/durandk/scratch_durandk/GenFAW600/VCF/202512_newvcf/Whole_genome_biallelic_max_5_missing_pruned.vcf.gz  \
# 	  --chr-set 29 \
# 	  --allow-extra-chr \
# 	  --pca \
# 	  --out Whole_genome_biallelic_max_5_missing_pruned_PCA



#############################################
# PCA autosomes
# #############################################

# plink2 --vcf /home/durandk/scratch_durandk/GenFAW600/VCF/202512_newvcf/Autosome_biallelic_max_5_missing_pruned.vcf.gz  \
# 	  --chr-set 29 \
# 	  --allow-extra-chr \
# 	  --pca \
# 	  --out Autosome_biallelic_max_5_missing_pruned_PCA
	 
#############################################
# PCA Z chromosome
# #############################################

plink2 --vcf /home/durandk/scratch_durandk/GenFAW600/VCF/202512_newvcf/Z_biallelic_max_5_missing_pruned.vcf.gz \
	  --chr-set 29 \
	  --allow-extra-chr \
	  --pca \
	  --out Z_biallelic_max_5_missing_pruned_PCA
