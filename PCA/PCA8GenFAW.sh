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
# PCA #--geno 0.2: eliminates SNPs with more than 20% missing, retaining those with at least 80% of data present.

# ==============================
# plink2 --bfile /storage/simple/users/durandk/scratch_durandk/GenFAW600/VCF/202512_newvcf/WholeGenome_biallelic_max80missing_pruned \
# 	  --chr-set 29 \
# 	  --allow-extra-chr \
# 	  --pca \
# 	  --out WholeGenome_biallelic_max80missing_pruned_PCA



#############################################
# PCA autosomes
#############################################

# plink2 --bfile /home/durandk/scratch_durandk/GenFAW600/VCF/202512_newvcf/Autosome_biallelic_max80missing_pruned \
# 	  --chr-set 29 \
# 	  --allow-extra-chr \
# 	  --pca \
# 	  --out Autosome_biallelic_max80missing_pruned_PCA
	 
#############################################
# PCA Z chromosome
#############################################

plink2 --bfile /home/durandk/scratch_durandk/GenFAW600/VCF/202512_newvcf/Z_biallelic_max80missing_pruned \
	  --chr-set 29 \
	  --allow-extra-chr \
	  --pca \
	  --out Z_biallelic_max80missing_pruned_PCA
