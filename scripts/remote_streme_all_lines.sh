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
#SBATCH --array=1-10%5

module load conda
eval "$(conda shell.bash hook)"
conda activate tf

# Define the path to the script
cd /home/l338m483/scratch/MEME_Test/All_lines

# Create array of FASTA files
FASTA_FILES=(*.fasta)
TOTAL_FILES=${#FASTA_FILES[@]}

# Check if array index is valid
if [ $SLURM_ARRAY_TASK_ID -gt $TOTAL_FILES ]; then
    echo "Array task ID $SLURM_ARRAY_TASK_ID exceeds number of files ($TOTAL_FILES)"
    exit 1
fi

# Get the current file (array index starts at 1, bash array starts at 0)
CURRENT_FILE=${FASTA_FILES[$((SLURM_ARRAY_TASK_ID-1))]}

echo "Processing file: $CURRENT_FILE"
echo "Task ID: $SLURM_ARRAY_TASK_ID"
echo "Job ID: $SLURM_ARRAY_JOB_ID"

# Run RepeatMasker for current file
echo "Running RepeatMasker on $CURRENT_FILE"
if RepeatMasker -species "Mimulus guttatus" "$CURRENT_FILE"; then
    echo "RepeatMasker completed successfully for $CURRENT_FILE"
else
    echo "RepeatMasker failed for $CURRENT_FILE"
    exit 1
fi

# Check if masked file was created
MASKED_FILE="${CURRENT_FILE}.masked"
if [ ! -f "$MASKED_FILE" ]; then
    echo "Masked file $MASKED_FILE not found"
    exit 1
fi

# Run STREME for current masked file
echo "Running STREME on $MASKED_FILE"
OUTPUT_DIR="streme_$(basename "$CURRENT_FILE" '.fasta')"

if streme --p "$MASKED_FILE" --nmotifs 200 --thresh 0.05 --minw 6 --maxw 20 -o "$OUTPUT_DIR"; then
    echo "STREME completed successfully for $MASKED_FILE"
    echo "Output saved to: $OUTPUT_DIR"
else
    echo "STREME failed for $MASKED_FILE"
    exit 1
fi

echo "Task completed successfully for $CURRENT_FILE"