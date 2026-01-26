#!/bin/bash
#SBATCH -p cpu-dedicated
#SBATCH --account=dedicated-cpu@dgimi-eha
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8
# ==============================================================================
# Script: filter_and_prune_GenFAW.sh
# Description: 
#   Filters a large VCF file to retain only high-quality biallelic SNPs with 
#   ≤5% missing genotypes, then performs linkage disequilibrium (LD) pruning 
#   using PLINK. The final output is a pruned, compressed, and indexed VCF.
#
# Input: 
#   Whole_genome_biallelic_max80missing.vcf.gz (already biallelic and filtered 
#   to ≤80% missingness)

# Tools required: PLINK (1.9 or 2.x), tabix
# ==============================================================================

# Chemin du VCF d'entrée
INPUT_VCF="/storage/simple/projects/faw_adaptation/Data_Backup/Merged_vcf/2025_GenFAW600/2025_GenFAW600/Whole_genome_biallelic_max80missing.vcf.gz"

# Préfixe pour les fichiers de sortie
PREFIX="GenFAW"


# ----------------------------- Step 1: bcftools filtering -----------------------------
source /home/durandk/miniconda3/etc/profile.d/conda.sh
conda activate bcftools
# # Create intermediate filtered VCF
# bcftools view "$INPUT_VCF" \
#     --threads 8 \
#     --include 'F_MISSING <= 0.05' \
#     --min-alleles 2 --max-alleles 2 \
#     --output-type z \
#     --output "${PREFIX}_max5miss_biallelic.vcf.gz"

# # Index the intermediate VCF (required for some tools and good practice)


# conda deactivate
# source /home/durandk/miniconda3/etc/profile.d/conda.sh
# conda activate bgzip_tabix

# tabix -p vcf "${PREFIX}_max5miss_biallelic.vcf.gz"
# ----------------------------- Step 2: PLINK LD Pruning -----------------------------
conda activate plink1.9

plink --vcf "${PREFIX}_max5miss_biallelic.vcf.gz" \
      --allow-extra-chr \
          --chr-set 29\
      --double-id \
      --indep-pairwise 50 10 0.2 \
      --out ${PREFIX}_prune

# Step 3: Extracting pruned variants and exporting to compressed VCF
plink --vcf "${PREFIX}_max5miss_biallelic.vcf.gz" \
      --allow-extra-chr \
            --double-id \
             --chr-set 29\
      --extract "${PREFIX}_prune.prune.in" \
      --recode vcf bgz \
      --out WholeGenome_biallelic_max5miss_pruned.vcf.gz

conda activate bgzip_tabix
# Step 4  index the final VCF
tabix -p vcf WholeGenome_biallelic_max5miss_pruned.vcf.gz


conda deactivate
source /home/durandk/miniconda3/etc/profile.d/conda.sh
conda activate bcftools
                 


#############################################
# 1. Extract autosomes (chromosomes 1 to 28)
#############################################

bcftools view \
    -r 1-28  WholeGenome_biallelic_max5miss_pruned.vcf.gz \
    -Oz -o Autosome_biallelic__max5miss_pruned.vcf.gz \
    --threads 8


#############################################
# 3. Extract chromosome Z (chromosome 29)
#############################################

bcftools view \
    -r 29  WholeGenome_biallelic_max5miss_pruned.vcf.gz \
    -Oz  -o Z_biallelic_max5miss_pruned.vcf.gz --threads 8

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
bcftools index Autosome_biallelic_max5miss_pruned.vcf.gz
bcftools index Z_biallelic_max5miss_pruned.vcf.gz
bcftools index -f merged.snps.max_5_missing_pruned.Treemix.vcf.gz









