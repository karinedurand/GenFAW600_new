#!/bin/bash
#SBATCH -p cpu-dedicated
#SBATCH --account=dedicated-cpu@dgimi-eha
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8

source /home/durandk/miniconda3/etc/profile.d/conda.sh
conda activate vcftools

#############################################
#Fix VCF header after PLINK pruning
#
# After LD pruning with PLINK, the VCF header contains a non-standard
# PLINK-specific line (##chrSet=...), which makes the file non-compliant
# with the VCF specification and prevents bcftools from parsing it.
## This step removes the problematic header line without modifying
# any variants or genotypes, and restores a valid VCF format.

#############################################
zgrep -v '^##chrSet=' WholeGenome_biallelic_max80missing_pruned.vcf.gz \
 | bgzip -c \
 > WholeGenome_biallelic_max80missing_pruned.fixed.vcf.gz
          
tabix -f -p vcf WholeGenome_biallelic_max80missing_pruned.fixed.vcf.gz

conda deactivate
source /home/durandk/miniconda3/etc/profile.d/conda.sh
conda activate bcftools
#############################################
# 1. Filter whole-genome VCF
#    - Keep only biallelic variants
#    - Keep variants with at least 95% present genotypes (F_MISSING < 0.05)
#    --keep variants with a maximum of 5% missing data (F_MISSING < 0.05)
#############################################

bcftools view   -i 'F_MISSING < 0.05'  WholeGenome_biallelic_max80missing_pruned.fixed.vcf.gz \
 -Oz  -o Whole_genome_biallelic_max_5_missing.pruned.vcf.gz \
 --threads 8                     


#############################################
# 2. Extract autosomes (chromosomes 1 to 28)
#############################################

bcftools view \
    -r 1-28  WholeGenome_biallelic_max80missing_pruned.fixed.vcf.gz \
    -Oz -o Autosome_biallelic_max_5_missing.pruned.vcf.gz \
    --threads 8


#############################################
# 3. Extract chromosome Z (chromosome 29)
#############################################

bcftools view \
    -r 29  WholeGenome_biallelic_max80missing_pruned.fixed.vcf.gz \
    -Oz  -o Z_biallelic_max_5_missing.pruned.vcf.gz --threads 8

#############################################
# 4. Filter whole-genome  with outgroup VCF
#############################################

bcftools view  -i 'F_MISSING < 0.05'  /storage/simple/projects/faw_adaptation/Data_Backup/Merged_vcf/2025_GenFAW600/2025_GenFAW600/merged.snps.filtered.renamed_chr1-29.miss0.8.vcf.gz  \
 -Oz  -o merged.snps.max_5_missing.Treemix.vcf.gz --threads 8              
zcat merged.snps.max_5_missing.Treemix.vcf.gz | grep -v "#" | wc -l 

#############################################
# 5. Index the resulting VCF files
#############################################

bcftools index Whole_genome_biallelic_max_5_missing.pruned.vcf.gz
bcftools index Autosome_biallelic_max_5_missing.pruned.vcf.gz
bcftools index Z_biallelic_max_5_missing.pruned.vcf.gz
bcftools index -f merged.snps.max_5_missing_pruned.Treemix.vcf.gz
