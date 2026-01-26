#!/bin/bash
#SBATCH -p cpu-dedicated
#SBATCH --account=dedicated-cpu@dgimi-eha
source /home/durandk/miniconda3/etc/profile.d/conda.sh
conda activate vcftools
#############################################
# /storage/simple/users/durandk/scratch_durandk/GenFAW600/VCF/202512_newvcf/merged.snps.max_5_missing.Treemix.vcf.gz
#  - Raw vcf rename chrom from 1-29
#    - Keep variants with at least 95% present genotypes (F_MISSING < 0.05) =  --keep variants with a maximum of 5% missing data (F_MISSING < 0.05)
############################################


sbatch vcf2treemix.sh /storage/simple/users/durandk/scratch_durandk/GenFAW600/VCF/202512_newvcf/merged.snps.max_5_missing.Treemix.vcf.gz pop_treemix.clust


