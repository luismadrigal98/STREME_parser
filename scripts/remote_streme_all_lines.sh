#!/bin/bash
#SBATCH --job-name=streme_array
#SBATCH --output=streme_array_%A_%a.output
#SBATCH --error=streme_array_%A_%a.error
#SBATCH --partition=eeb,kelly,kucg
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=4g
#SBATCH --time=10-00:00:00
#SBATCH --mail-user=madrigalrocalj@ku.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --array=1-2

# Genome-generic STREME preparation as a SLURM array job.
# Each array task prepares one genome from a manifest:
#   genome -> promoters -> mask -> background -> STREME (streme_<genome>/).
#
# Submit AFTER editing the variables below:
#   sbatch --array=1-<N> scripts/remote_streme_all_lines.sh
# where <N> is the number of genome rows in the manifest.

set -euo pipefail

module load conda
eval "$(conda shell.bash hook)"
conda activate tf

# ----- Edit these -----
REPO_DIR="${REPO_DIR:-$HOME/MEME_related}"          # path to this repository
MANIFEST="${MANIFEST:-genomes.tsv}"                  # genome  fasta  annotation [expression]
OUTPUT_ROOT="${OUTPUT_ROOT:-prepared}"               # where streme_<genome>/ dirs are written
UPSTREAM="${UPSTREAM:-1000}"
MASKER="${MASKER:-repeatmasker}"                     # repeatmasker | dust | none
SPECIES="${SPECIES:-}"                               # RepeatMasker species (optional)
THREADS="${SLURM_CPUS_PER_TASK:-10}"                 # threads for the heavy steps
# ----------------------

# Pick the manifest row for this array task (skip the header line).
ROW_NUM=$((SLURM_ARRAY_TASK_ID + 1))
LINE=$(sed -n "${ROW_NUM}p" "$MANIFEST")
if [ -z "$LINE" ]; then
    echo "No manifest row for array task $SLURM_ARRAY_TASK_ID (line $ROW_NUM of $MANIFEST)"
    exit 1
fi

GENOME=$(echo "$LINE" | cut -f1)
FASTA=$(echo "$LINE" | cut -f2)
ANNOT=$(echo "$LINE" | cut -f3)

echo "Task $SLURM_ARRAY_TASK_ID -> genome=$GENOME fasta=$FASTA annotation=$ANNOT"

SPECIES_ARG=()
if [ -n "$SPECIES" ]; then
    SPECIES_ARG=(--species "$SPECIES")
fi

python "$REPO_DIR/pipelines/streme_pipeline.py" prepare \
    --genome "$GENOME" --fasta "$FASTA" --annotation "$ANNOT" \
    --upstream "$UPSTREAM" --mask "$MASKER" "${SPECIES_ARG[@]}" \
    --threads "$THREADS" --output "$OUTPUT_ROOT"

echo "Done: $GENOME -> $OUTPUT_ROOT/streme_$GENOME"
echo "After all tasks finish, consolidate with:"
echo "  python $REPO_DIR/cli_tools/streme_sites_consolidator.py consolidate $OUTPUT_ROOT --output outputs/consolidated_streme_sites"
