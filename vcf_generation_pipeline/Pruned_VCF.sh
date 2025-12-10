#!/bin/bash
#SBATCH -p cpu-dedicated
#SBATCH --account=dedicated-cpu@dgimi-eha
#SBATCH --mem=40G


source /home/durandk/miniconda3/etc/profile.d/conda.sh
conda activate plink2

plink2 \
    --vcf Whole_genome_biallelic_max80missing.vcf.gz \
    -set-all-var-ids @:#:$r:$a \
    --rm-dup exclude-all \
    --make-bed \
    --chr-set 29 \
    --allow-extra-chr \
    --out WholeGenome_filtered


#############################################
# 4. LD pruning
#############################################

plink2 \
    --bfile WholeGenome_filtered \
    --chr-set 29 \
    --allow-extra-chr \
    --indep-pairwise 50 5 0.2 \
    --out LDpruned


#############################################
# 5. Extract the pruned SNP set
#############################################

plink2 \
    --bfile WholeGenome_filtered \
    --chr-set 29 \
    --allow-extra-chr \
    --extract LDpruned.prune.in \
    --make-bed \
    --out WholeGenome_biallelic_max80missing_pruned 

plink2 --bfile WholeGenome_biallelic_max80missing_pruned --chr-set 29 --allow-extra-chr  --recode vcf --out WholeGenome_biallelic_max80missing_pruned

conda activate vcftools

bgzip .c  WholeGenome_biallelic_max80missing_pruned.vcf > WholeGenome_biallelic_max80missing_pruned.vcf.gz
tabix -p vcf WholeGenome_biallelic_max80missing_pruned.vcf.gz
