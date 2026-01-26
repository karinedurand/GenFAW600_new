#!/bin/bash
#SBATCH -p cpu-dedicated
#SBATCH --account=dedicated-cpu@dgimi-eha
#SBATCH --mem 50G

source /home/durandk/miniconda3/etc/profile.d/conda.sh
conda activate vcftools
VCF="/storage/simple/projects/faw_adaptation/Data_Backup/Merged_vcf/2025_GenFAW600/2025_GenFAW600/VCF_notpruned/GenFAW_max5miss_biallelic.vcf.gz"

vcftools --gzvcf $VCF  --weir-fst-pop Inv_east.pop  --weir-fst-pop Inv_west.pop --fst-window-size 500000 --fst-window-step 50000  --out Inv_east_Inv_west
 
vcftools --gzvcf $VCF  --weir-fst-pop Inv_east.pop  --weir-fst-pop Mex.pop --fst-window-size 500000 --fst-window-step 50000  --out Inv_east_Mex
 
vcftools --gzvcf $VCF  --weir-fst-pop Inv_east.pop --weir-fst-pop Nat_C.pop --fst-window-size 500000 --fst-window-step 50000 --out Inv_east_Nat_C
 
vcftools --gzvcf $VCF  --weir-fst-pop Inv_east.pop --weir-fst-pop Nat_R.pop --fst-window-size 500000 --fst-window-step 50000 --out Inv_east_Nat_R

vcftools --gzvcf $VCF  --weir-fst-pop Zam.pop --weir-fst-pop Inv_east.pop --fst-window-size 500000 --fst-window-step 50000 --out Zam_Inv_east

vcftools --gzvcf $VCF  --weir-fst-pop Mex.pop --weir-fst-pop Inv_west.pop --fst-window-size 500000 --fst-window-step 50000 --out Mex_Inv_west

vcftools --gzvcf $VCF  --weir-fst-pop Inv_west.pop --weir-fst-pop Nat_C.pop --fst-window-size 500000 --fst-window-step 50000 --out Inv_west_Nat_C

vcftools --gzvcf $VCF  --weir-fst-pop Inv_west.pop --weir-fst-pop Nat_R.pop --fst-window-size 500000 --fst-window-step 50000 --out Inv_west_Nat_R

vcftools --gzvcf $VCF  --weir-fst-pop Zam.pop --weir-fst-pop Inv_west.pop --fst-window-size 500000 --fst-window-step 50000 --out Zam_Inv_west

vcftools --gzvcf $VCF  --weir-fst-pop Mex.pop --weir-fst-pop Nat_C.pop --fst-window-size 500000 --fst-window-step 50000 --out Mex_Nat_C

vcftools --gzvcf $VCF  --weir-fst-pop Mex.pop --weir-fst-pop Nat_R.pop --fst-window-size 500000 --fst-window-step 50000 --out Mex_Nat_R

vcftools --gzvcf $VCF  --weir-fst-pop Zam.pop  --weir-fst-pop Mex.pop --fst-window-size 500000 --fst-window-step 50000  --out Zam_Mex
  
vcftools --gzvcf $VCF  --weir-fst-pop Nat_C.pop --weir-fst-pop Nat_R.pop --fst-window-size 500000 --fst-window-step 50000 --out Nat_C_Nat_R

vcftools --gzvcf $VCF  --weir-fst-pop Zam.pop --weir-fst-pop Nat_C.pop --fst-window-size 500000 --fst-window-step 50000 --out Zam_Nat_C

vcftools --gzvcf $VCF  --weir-fst-pop Zam.pop --weir-fst-pop Nat_R.pop --fst-window-size 500000 --fst-window-step 50000 --out Zam_Nat_R


