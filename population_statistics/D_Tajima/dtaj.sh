#!/bin/bash
#SBATCH -p cpu-dedicated
#SBATCH --account=dedicated-cpu@dgimi-eha
#SBATCH --array=6
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00

source /home/durandk/miniconda3/etc/profile.d/conda.sh
conda activate bgzip_tabix
cp /storage/simple/projects/faw_adaptation/Data_Backup/Merged_vcf/2025_GenFAW600/2025_GenFAW600/VCF_notpruned/GenFAW_max5miss_biallelic.vcf.gz*   .
# --- Input VCF ---

VCF="GenFAW_max5miss_biallelic.vcf.gz"
# --- Population files ---
POP=(
"Inv_west.pop"
"Mex.pop"
"Nat_C.pop"
"Nat_R.pop"
"Zam.pop"
"Inv_east.pop"
)

# --- Output directory ---
OUTDIR="tajimaD_results_new"
mkdir -p "$OUTDIR"

# --- Environment ---
source /home/durandk/miniconda3/etc/profile.d/conda.sh
conda activate vcftools

# --- Get current population ---
CURRENT_POP="${POP[$((SLURM_ARRAY_TASK_ID - 1))]}"
POP_NAME=$(basename "$CURRENT_POP" .pop)  
OUTPREFIX="${OUTDIR}/${POP_NAME}_TajimaD"


vcftools --gzvcf "$VCF" \
  --keep "$CURRENT_POP" \
  --TajimaD 100000 \
  --out "$OUTPREFIX"


