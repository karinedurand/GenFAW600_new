#!/bin/bash
#SBATCH -p cpu-dedicated
#SBATCH --account=dedicated-cpu@dgimi-eha
#SBATCH --mem=80G
#SBATCH -c 1


source /home/durandk/miniconda3/etc/profile.d/conda.sh
conda activate vcftools

VCF="/storage/simple/projects/faw_adaptation/Data_Backup/Merged_vcf/2025_GenFAW600/2025_GenFAW600/VCF_notpruned/GenFAW_max5miss_biallelic.vcf.gz"
/
for i in {1..29}
do
  vcftools --gzvcf $VCF --LROH --chr $i --out ${i}_roh_result
done

conda deactivate
