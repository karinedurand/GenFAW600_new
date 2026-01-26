#!/bin/bash
#SBATCH -p cpu-dedicated
#SBATCH --account=dedicated-cpu@dgimi-eha
#SBATCH --mem=16G
#SBATCH --cpus-per-task=4
#SBATCH --job-name=treemix_optm
#SBATCH --array=0-6        
#SBATCH --time=24:00:00    


module load bioinfo-ifb
module load treemix/1.13

# ====================== PARAMÈTRES ======================
NITER=10                   
M=${SLURM_ARRAY_TASK_ID}    
INPUT="merged.snps.max_5_missing.Treemix.popsNMF.frq.gz"
OUTDIR="treemix_optm_runs" 
SEED_BASE=123456

# Création du dossier de sortie
mkdir -p ${OUTDIR}

echo "=== Lancement TreeMix pour m=${M} avec ${NITER} replicates ==="

# Boucle sur les replicates
for i in $(seq 1 ${NITER}); do
    # Seed unique par replicate et par m → variance garantie
    SEED=$((SEED_BASE + i * 100 + M * 1000))

    # Légère variation de -k pour plus de robustesse et variance
    K=$((500 + (i - 1) * 50))  # 500, 550, 600, ..., 950

    # Nom de sortie STANDARD pour OptM
    OUT="${OUTDIR}/treemix.m${M}.rep${i}"

    echo "Run ${i}/${NITER} | m=${M} | rep=${i} | seed=${SEED} | k=${K} | output=${OUT}"

    treemix \
        -i ${INPUT} \
        -o ${OUT} \
        -root slitura \
        -global \
        -k ${K} \
        -seed ${SEED} \
        -m ${M}

done

echo "=== Terminé pour m=${M} ==="
