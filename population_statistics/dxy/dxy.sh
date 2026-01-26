#!/bin/bash
#SBATCH -p cpu-dedicated
#SBATCH --account=dedicated-cpu@dgimi-eha
#SBATCH --array=1-1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G

cd /storage/simple/users/durandk/scratch_durandk/GenFAW600/dxy


# Définir les comparaisons de populations
POP_PAIRS=(
  "Inv_east_Inv_west"
  "Inv_east_Mex"
  "Inv_east_Nat_C"
  "Inv_east_Nat_R"
  "Inv_east_Zam"
  "Inv_west_Mex"
 "Inv_west_Nat_C"
  "Inv_west_Nat_R"
  "Inv_west_Zam"
  "Mex_Nat_C"
  "Mex_Nat_R"
  "Mex_Zam"
 "Nat_C_Nat_R"
  "Nat_C_Zam"
  "Nat_R_Zam"
)

# Récupération du chromosome à traiter à partir de l'index d'array
CHR=${SLURM_ARRAY_TASK_ID}

# Fichier VCF source
VCF_SRC="/storage/simple/projects/faw_adaptation/Data_Backup/Merged_vcf/2025_GenFAW600/2025_GenFAW600/VCF_notpruned/GenFAW_max5miss_biallelic.vcf.gz"

# Extraction du chromosome correspondant
source /home/durandk/miniconda3/etc/profile.d/conda.sh
conda activate bgzip_tabix

VCF_TMP=${CHR}.vcf.gz
tabix -h "$VCF_SRC" "$CHR" | gzip -c > "$VCF_TMP"

conda deactivate
module load bioinfo-cirad
module load anaconda/python3.8
# Calcul du Dxy pour chaque paire de populations
for PAIR in "${POP_PAIRS[@]}"; do
    POP1=${PAIR%%_*}
    POP2=${PAIR##*_}
    OUTNAME="${POP1}_${POP2}_chr${CHR}.dxy"
   # python3 Dxy_calculate -v "$VCF_TMP" -p "$PAIR" -o "$OUTNAME" -w 100000
   python3 Dxy_calculate -v   ${CHR}.vcf.gz -p "$PAIR" -o "$OUTNAME" -w 100000 -s 10000
done






