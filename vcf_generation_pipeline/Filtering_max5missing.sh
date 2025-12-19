#!/bin/bash
#SBATCH -p cpu-dedicated
#SBATCH --account=dedicated-cpu@dgimi-eha
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8


source /home/durandk/miniconda3/etc/profile.d/conda.sh
conda activate bcftools
#############################################
# 1. Filter whole-genome VCF
#    - Keep only biallelic variants
#    - Keep variants with at least 95% present genotypes (F_MISSING < 0.05)
#    --keep variants with a maximum of 5% missing data (F_MISSING < 0.05)
#############################################

bcftools view     -m2 -M2 -i 'F_MISSING < 0.05' GenFAW.snps.renamed_chr1-29.vcf.gz \
 -Oz  -o Whole_genome_biallelic_max_5_missing.vcf.gz \
 --threads 8                     


#############################################
# 2. Extract autosomes (chromosomes 1 to 28)
#############################################

bcftools view \
    -r 1-28 Whole_genome_biallelic_max_5_missing.vcf.gz  \
    -Oz -o Autosome_biallelic_max_5_missing.vcf.gz \
    --threads 8


#############################################
# 3. Extract chromosome Z (chromosome 29)
#############################################

bcftools view \
    -r 29 Whole_genome_biallelic_max_5_missing.vcf.gz  \
    -Oz  -o Z_biallelic_max_5_missing.vcf.gz --threads 8


#############################################
# 4. Index the resulting VCF files
#############################################

bcftools index Whole_genome_biallelic_max_5_missing.vcf.gz 
bcftools index Autosome_biallelic_max_5_missing.vcf.gz
bcftools index Z_biallelic_max_5_missing.vcf.gz
