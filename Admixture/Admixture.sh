#!/bin/bash
#SBATCH -p cpu-dedicated
#SBATCH --account=dedicated-cpu@dgimi-eha
#SBATCH --nodes=1  
#SBATCH --cpus-per-task=8 
#SBATCH --mem=120G
#SBATCH -J admixture

# ------------------------------------------------------------
# PLINK filtering (geno only, no MAF filter)  - Remove SNPs with more than 5% missing genotypes (--geno 0.05)
# ------------------------------------------------------------

plink \
  --bfile /home/durandk/scratch_durandk/GenFAW600/VCF/202512_newvcf/WholeGenome_biallelic_max80missing_pruned \
  -chr-set 29 \
  --allow-extra-chr \
  --geno 0.05 \
  --make-bed \
  --out admix_geno0.05


# ------------------------------------------------------------
# 2. ADMIXTURE analyses
# ------------------------------------------------------------

K=${SLURM_ARRAY_TASK_ID}

# /storage/simple/projects/faw_adaptation/programs/scripts_SNP/admixture/admixture \
#   -j8 \
#   --cv \
#   admix_geno0.05.bed \
#   ${K} \
#   > log_admixture_geno0.05_K${K}.log 2>&1

