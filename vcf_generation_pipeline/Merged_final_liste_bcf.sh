#!/bin/bash
#SBATCH -p cpu-dedicated
#SBATCH --account=dedicated-cpu@dgimi-eha
#SBATCH --mem=40G
#SBATCH --cpus-per-task=8


source /home/durandk/miniconda3/etc/profile.d/conda.sh
conda activate bcftools

#-------------------------------
# 1. Variables
# -------------------------------
LIST="liste_bcf_amerged.txt"   # ta liste d’échantillons
DIR="/home/durandk/scratch_durandk/all_bcf"   # répertoire des liens BCF



# 3. Merge BCF 
# -------------------------------


bcftools merge \
   --file-list "$LIST" \
   -Oz \
   -o  GenFAW.vcf.gz \
   --threads 4

# -------------------------------
# 4. SNPs only
# -------------------------------

 bcftools view \
     -v snps \
     GenFAW.vcf.gz  \
     -Oz -o GenFAW.snps.vcf.gz \
      --threads 8

# -------------------------------
# 5. Index
# -------------------------------

bcftools index GenFAW.snps.vcf.gz

# -------------------------------
# 6. Rename chromosomes
# -------------------------------
 bcftools annotate --rename-chrs rename.txt GenFAW.snps.vcf.gz \
   -Oz -o GenFAW.snps.vcf.renamed.vcf.gz

bcftools index GenFAW.snps.vcf.renamed.vcf.gz
bcftools view -r 1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29 \
   GenFAW.snps.vcf.renamed.vcf.gz -Oz \
   -o GenFAW.snps.renamed_chr1-29.vcf.gz --threads 4
bcftools index GenFAW.snps.renamed_chr1-29.vcf.gz




