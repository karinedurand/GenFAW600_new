#!/bin/bash
#SBATCH -p cpu-dedicated
#SBATCH --account=dedicated-cpu@dgimi-eha
#SBATCH --mem=40G

#############################################
# 1. Prepare whole VCF for LD pruning
#############################################
source /home/durandk/miniconda3/etc/profile.d/conda.sh
conda activate plink2

# plink2 \
#     --vcf Whole_genome_biallelic_max80missing.vcf.gz \
#     -set-all-var-ids @:#:$r:$a \
#     --rm-dup exclude-all \
#     --make-bed \
#     --chr-set 29 \
#     --allow-extra-chr \
#     --out WholeGenome_filtered


# #############################################
# # 2. LD pruning 
# #############################################

# plink2 \
#     --bfile WholeGenome_filtered \
#     --chr-set 29 \
#     --allow-extra-chr \
#     --indep-pairwise 50 5 0.2 \
#     --out LDpruned


# #############################################
# # 3. Extract the pruned SNP set
# #############################################

# plink2 \
#     --bfile WholeGenome_filtered \
#     --chr-set 29 \
#     --allow-extra-chr \
#     --extract LDpruned.prune.in \
#     --make-bed \
#     --out 

# plink2 --bfile WholeGenome_biallelic_max80missing_pruned --chr-set 29 --allow-extra-chr  --recode vcf --out WholeGenome_biallelic_max80missing_pruned

# conda activate vcftools

# bgzip -c  WholeGenome_biallelic_max80missing_pruned.vcf > 
# tabix -p vcf WholeGenome_biallelic_max80missing_pruned.vcf.gz
# conda deactivate

#############################################
# 1. Prepare Autosome VCF for LD pruning
#############################################
# source /home/durandk/miniconda3/etc/profile.d/conda.sh
# conda activate plink2
# plink2 \
#     --vcf Autosome_biallelic_max80missing.vcf.gz \
#     -set-all-var-ids @:#:$r:$a \
#     --rm-dup exclude-all \
#     --make-bed \
#     --chr-set 29 \
#     --allow-extra-chr \
#     --out Autosome_filtered


# #############################################
# # 4. LD pruning
# #############################################

# plink2 \
#     --bfile Autosome_filtered \
#     --chr-set 29 \
#     --allow-extra-chr \
#     --indep-pairwise 50 5 0.2 \
#     --out LDpruned


# #############################################
# # 5. Extract the pruned SNP set
# #############################################

# plink2 \
#     --bfile Autosome_filtered \
#     --chr-set 29 \
#     --allow-extra-chr \
#     --extract LDpruned.prune.in \
#     --make-bed \
#     --out WholeGenome_biallelic_max80missing_pruned 

# plink2 --bfile WholeGenome_biallelic_max80missing_pruned --chr-set 29 --allow-extra-chr  --recode vcf --out WholeGenome_biallelic_max80missing_pruned

# conda activate vcftools

# bgzip -c  WholeGenome_biallelic_max80missing_pruned.vcf > WholeGenome_biallelic_max80missing_pruned.vcf.gz
# tabix -p vcf WholeGenome_biallelic_max80missing_pruned.vcf.gz

#############################################
# 1. Prepare Z VCF for LD pruning
#############################################
source /home/durandk/miniconda3/etc/profile.d/conda.sh
conda activate plink2
plink2 \
    --vcf Z_biallelic_max80missing.vcf.gz \
    -set-all-var-ids @:#:$r:$a \
    --rm-dup exclude-all \
    --make-bed \
    --chr-set 29 \
    --allow-extra-chr \
    --out Z_filtered


#############################################
# 4. LD pruning
#############################################

plink2 \
    --bfile Z_filtered\
    --chr-set 29 \
    --allow-extra-chr \
    --indep-pairwise 50 5 0.2 \
    --out LDpruned


#############################################
# 5. Extract the pruned SNP set
#############################################

plink2 \
    --bfile Z_filtered \
    --chr-set 29 \
    --allow-extra-chr \
    --extract LDpruned.prune.in \
    --make-bed \
    --out Z_biallelic_max80missing_pruned 

plink2 --bfile Z_biallelic_max80missing_pruned --chr-set 29 --allow-extra-chr  --recode vcf --out Z_biallelic_max80missing_pruned

conda activate vcftools

bgzip -c  Z_biallelic_max80missing_pruned.vcf > Z_biallelic_max80missing_pruned.vcf.gz
tabix -p vcf Z_biallelic_max80missing_pruned.vcf.gz
cp Z_biallelic_max80missing_pruned.vcf.gz* /storage/simple/projects/faw_adaptation/Data_Backup/Merged_vcf/2025_GenFAW600/2025_GenFAW600/

