#!/bin/bash
#SBATCH --job-name=prep_P_eatonii
#SBATCH --cpus-per-task=16
#SBATCH --mem=96g
#SBATCH --time=2-00:00:00
#SBATCH --partition=eeb,kelly,kucg
#SBATCH --output=prep_P_eatonii_%j.out
#SBATCH --error=prep_P_eatonii_%j.err
#SBATCH --mail-user=madrigalrocalj@ku.edu
#SBATCH --mail-type=END,FAIL

# Stages 2+3 of Pipeline_P_eatonii.yaml:
#  - resolve the 1.x library to a stable symlink (repeats/P_eatonii/library.fa)
#  - extract promoters → RepeatMasker (custom lib) → TRF → order-3 background
#    → STREME (50 motifs) → FIMO (q ≤ 0.05, max-strand) inline
# STREME is the heavy step memory-wise (~30-50 GB peak on ~27k 2-kb promoters).

module load conda
eval "$(conda shell.bash hook)"
conda activate PyR

cd /kuhpc/scratch/kelly/l338m483/MEME_Penstemon

# Stage 2: resolve the 1.x repeat library to a stable path the rest of the
# pipeline can reference without depending on the timestamped RM_<date>/ dir.
LIB=$(ls -t repeats/P_eatonii/RM_*/consensi.fa.classified repeats/P_eatonii/RM_*/consensi.fa 2>/dev/null | head -1)
if [ -z "${LIB:-}" ] || [ ! -s "$LIB" ]; then
    echo "ERROR: no RepeatModeler library found under repeats/P_eatonii/RM_*/" >&2
    echo "       Did the RepeatModeler job finish successfully?" >&2
    exit 1
fi
ln -sfn "$(realpath "$LIB")" repeats/P_eatonii/library.fa
echo "library -> $(realpath repeats/P_eatonii/library.fa)"

# Stage 3: prepare + STREME + FIMO inline.
sh /home/l338m483/bin/STREME_parser/bin/streme-parser prepare \
    --genome P_eatonii \
    --fasta penstemon_eatonii/PeChr.BYU.final.fa \
    --annotation penstemon_eatonii/penstemon_eatonii_final_nuclear_cp.gff3 \
    --contig-pattern '^PeChr' \
    --upstream 2000 \
    --avoid-overlap \
    --mask dust \
    --extra-mask trf \
    --background-order 3 \
    --nmotifs 50 \
    --threads "$SLURM_CPUS_PER_TASK" \
    --output prepared_penstemon_eatonii \
    --run-fimo --fimo-qv-thresh --fimo-thresh 0.05 --fimo-max-strand

# Diagnostic: how much got masked? Target for plant promoters: >= 30%.
awk '!/^>/{n+=length($0); for(i=1;i<=length($0);i++){c=substr($0,i,1); if(c=="N"||c~/[a-z]/) m++}} END{printf "masked fraction: %.3f\n", m/n}' \
    prepared_penstemon_eatonii/P_eatonii_prep/P_eatonii_promoters.masked.trf
