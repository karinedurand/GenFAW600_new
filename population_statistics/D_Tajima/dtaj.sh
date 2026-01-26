#!/bin/bash
#SBATCH -p cpu-dedicated
#SBATCH --account=dedicated-cpu@dgimi-eha
#SBATCH --array=6
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --time=48:00:00

source /home/durandk/miniconda3/etc/profile.d/conda.sh
conda activate bgzip_tabix
#cp /home/durandk/scratch_durandk/GenFAW600/VCF/202512_newvcf/Whole_genome_biallelic_max_5_missing_pruned.vcf.gz*   .

# zcat Whole_genome_biallelic_max_5_missing_pruned.vcf.gz\
#   | sed 's/^##fileformat=VCFv4.3/##fileformat=VCFv4.2/' \
#    |bgzip -c > Whole_genome_biallelic_max_5_missing_pruned_v4.2.vcf.gz


# tabix -p vcf Whole_genome_biallelic_max_5_missing_pruned_v4.2.vcf.gz

# --- Input VCF ---

VCF="Whole_genome_biallelic_max_5_missing_pruned_v4.2.vcf.gz"
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


