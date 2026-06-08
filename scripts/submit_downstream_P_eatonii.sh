#!/bin/bash
#SBATCH --job-name=down_P_eatonii
#SBATCH --cpus-per-task=4
#SBATCH --mem=16g
#SBATCH --time=04:00:00
#SBATCH --partition=eeb,kelly,kucg
#SBATCH --output=down_P_eatonii_%j.out
#SBATCH --error=down_P_eatonii_%j.err
#SBATCH --mail-user=madrigalrocalj@ku.edu
#SBATCH --mail-type=END,FAIL

# Stages 4-6 of Pipeline_P_eatonii.yaml:
#  4. TOMTOM: STREME motifs vs A. thaliana TF PWMs (PlantTFDB)
#  5. consolidate + validate the within-genome motif clustering
#  6. exploratory motif co-occurrence network + gene clustering
# All quick; chained in one short job.

module load conda
eval "$(conda shell.bash hook)"
conda activate PyR

cd /kuhpc/scratch/kelly/l338m483/MEME_Penstemon

# Stage 4: TOMTOM
streme-parser annotate prepared_penstemon_eatonii/ \
    --target-db databases/Ath_TF_binding_motifs.meme \
    --thresh 0.1 --jobs 1

# Stage 5: consolidate + validate
streme-parser consolidate prepared_penstemon_eatonii/ \
    --output outputs/consolidated_P_eatonii \
    --threshold 0.75

streme-parser validate outputs/consolidated_P_eatonii.tsv

# Stage 6: motif co-occurrence network + gene clustering
streme-parser network outputs/consolidated_P_eatonii.tsv \
    --output network_P_eatonii/ \
    --min-motif-sites 20 \
    --min-jaccard 0.15 --min-lift 2.5 --fdr 0.01 \
    --n-clusters 30

echo
echo "All downstream stages complete."
echo "  TOMTOM   -> prepared_penstemon_eatonii/tomtom_P_eatonii/tomtom.tsv"
echo "  Consol.  -> outputs/consolidated_P_eatonii.tsv"
echo "  Network  -> network_P_eatonii/motif_cooccurrence_network.graphml"
