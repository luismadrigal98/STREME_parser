#!/bin/bash
#SBATCH --job-name=rm_P_eatonii
#SBATCH --cpus-per-task=16
#SBATCH --mem=32g
#SBATCH --time=5-00:00:00
#SBATCH --partition=eeb,kelly,kucg
#SBATCH --output=rm_P_eatonii_%j.out
#SBATCH --error=rm_P_eatonii_%j.err
#SBATCH --mail-user=madrigalrocalj@ku.edu
#SBATCH --mail-type=END,FAIL

# Stage 1 of Pipeline_P_eatonii.yaml: build the custom RepeatMasker library
# with RepeatModeler 1.x (--legacy in the wrapper switches to -pa, -engine ncbi,
# drops -LTRStruct). The library lands at RM_<timestamp>/consensi.fa.classified;
# Stage 2 (in scripts/submit_prepare_P_eatonii.sh) symlinks it to a stable path.

module load conda
eval "$(conda shell.bash hook)"
conda activate PyR

cd /kuhpc/scratch/kelly/l338m483/MEME_Penstemon

python /home/l338m483/bin/STREME_parser/cli_tools/genome_prep.py model-repeats \
    penstemon_eatonii/PeChr.BYU.final.fa \
    --output-dir repeats/P_eatonii \
    --name P_eatonii \
    --threads "$SLURM_CPUS_PER_TASK" \
    --legacy

echo "Done. Library at: $(ls -t repeats/P_eatonii/RM_*/consensi.fa.classified | head -1)"
