#!/bin/bash
# Consolidate STREME results from multiple line analyses
# This script processes the output from your array job

set -e  # Exit on any error

echo "=== STREME Results Consolidation Script ==="
echo "Processing results from individual line analyses..."

# Define directories
WORK_DIR="/home/l338m483/scratch/MEME_Test/All_lines"
CONSOLIDATION_DIR="$WORK_DIR/consolidated_results"

# Create consolidation directory
mkdir -p "$CONSOLIDATION_DIR"

cd "$WORK_DIR"

# Find all STREME output directories
STREME_DIRS=(streme_*)
NUM_DIRS=${#STREME_DIRS[@]}

echo "Found $NUM_DIRS STREME output directories"

if [ $NUM_DIRS -eq 0 ]; then
    echo "No STREME output directories found!"
    exit 1
fi

# List all directories for verification
echo "STREME directories found:"
for dir in "${STREME_DIRS[@]}"; do
    echo "  $dir"
done

# Create summary file
SUMMARY_FILE="$CONSOLIDATION_DIR/analysis_summary.txt"
echo "Creating analysis summary: $SUMMARY_FILE"

cat > "$SUMMARY_FILE" << EOF
# STREME Analysis Summary - All Lines
# Generated on $(date)
# Total lines analyzed: $NUM_DIRS
# 
# Directory Structure:
EOF

# Process each STREME directory
echo ""
echo "Processing individual STREME results..."

for i in "${!STREME_DIRS[@]}"; do
    dir="${STREME_DIRS[i]}"
    line_num=$((i + 1))
    
    echo "Processing $dir (Line $line_num)..."
    
    # Add to summary
    echo "Line_$line_num: $dir" >> "$SUMMARY_FILE"
    
    # Check if streme.txt exists
    if [ -f "$dir/streme.txt" ]; then
        echo "  ✓ streme.txt found"
        
        # Extract motif count
        motif_count=$(grep -c "^MOTIF" "$dir/streme.txt" 2>/dev/null || echo "0")
        echo "  ✓ Found $motif_count motifs"
        
        # Copy important files with line prefix
        cp "$dir/streme.txt" "$CONSOLIDATION_DIR/line_${line_num}_streme.txt"
        
        if [ -f "$dir/streme.html" ]; then
            cp "$dir/streme.html" "$CONSOLIDATION_DIR/line_${line_num}_streme.html"
        fi
        
    else
        echo "  ✗ streme.txt not found in $dir"
        echo "Line_$line_num: $dir - MISSING streme.txt" >> "$SUMMARY_FILE"
    fi
done

echo ""
echo "=== Pair-wise Consolidation ==="

# Create pair-wise analysis
PAIRS_DIR="$CONSOLIDATION_DIR/pair_analysis"
mkdir -p "$PAIRS_DIR"

echo "Creating pair-wise consolidation..."

# Process lines in pairs
for ((i=0; i<NUM_DIRS; i+=2)); do
    pair_num=$(((i/2) + 1))
    line1=$((i + 1))
    line2=$((i + 2))
    
    echo "Processing Pair $pair_num: Line $line1 and Line $line2"
    
    pair_file="$PAIRS_DIR/pair_${pair_num}_lines_${line1}_${line2}.txt"
    
    cat > "$pair_file" << EOF
# Pair Analysis: Lines $line1 and $line2
# Generated on $(date)
#
# ==============================================

## LINE $line1 ANALYSIS
EOF
    
    # Add Line 1 results if available
    if [ -f "$CONSOLIDATION_DIR/line_${line1}_streme.txt" ]; then
        echo "Source: ${STREME_DIRS[$i]}" >> "$pair_file"
        echo "" >> "$pair_file"
        
        # Extract motif information
        grep -A 5 "^MOTIF" "$CONSOLIDATION_DIR/line_${line1}_streme.txt" >> "$pair_file" 2>/dev/null || echo "No motifs found" >> "$pair_file"
    else
        echo "No results available for Line $line1" >> "$pair_file"
    fi
    
    echo "" >> "$pair_file"
    echo "## LINE $line2 ANALYSIS" >> "$pair_file"
    
    # Add Line 2 results if available and line2 exists
    if [ $line2 -le $NUM_DIRS ] && [ -f "$CONSOLIDATION_DIR/line_${line2}_streme.txt" ]; then
        echo "Source: ${STREME_DIRS[$((i+1))]}" >> "$pair_file"
        echo "" >> "$pair_file"
        
        # Extract motif information
        grep -A 5 "^MOTIF" "$CONSOLIDATION_DIR/line_${line2}_streme.txt" >> "$pair_file" 2>/dev/null || echo "No motifs found" >> "$pair_file"
    else
        echo "No results available for Line $line2 (single line or missing)" >> "$pair_file"
    fi
    
    echo "" >> "$pair_file"
    echo "## PAIR COMPARISON" >> "$pair_file"
    echo "# TODO: Add motif similarity analysis here" >> "$pair_file"
    echo "# Consider using TOMTOM for motif comparison" >> "$pair_file"
    
done

# Create master consolidation file
MASTER_FILE="$CONSOLIDATION_DIR/master_consolidation.txt"
echo "Creating master consolidation file: $MASTER_FILE"

cat > "$MASTER_FILE" << EOF
# Master STREME Results Consolidation
# Generated on $(date)
# Analysis of $NUM_DIRS lines processed individually
#
# ==============================================

## SUMMARY STATISTICS
Total lines analyzed: $NUM_DIRS
Total pairs created: $(((NUM_DIRS + 1) / 2))

## MOTIF COUNTS PER LINE
EOF

# Add motif counts
for i in "${!STREME_DIRS[@]}"; do
    line_num=$((i + 1))
    if [ -f "$CONSOLIDATION_DIR/line_${line_num}_streme.txt" ]; then
        motif_count=$(grep -c "^MOTIF" "$CONSOLIDATION_DIR/line_${line_num}_streme.txt" 2>/dev/null || echo "0")
        echo "Line $line_num: $motif_count motifs" >> "$MASTER_FILE"
    else
        echo "Line $line_num: No results" >> "$MASTER_FILE"
    fi
done

echo "" >> "$MASTER_FILE"
echo "## FILES GENERATED" >> "$MASTER_FILE"
echo "Individual results: $CONSOLIDATION_DIR/line_*_streme.txt" >> "$MASTER_FILE"
echo "Pair analyses: $PAIRS_DIR/pair_*_lines_*.txt" >> "$MASTER_FILE"
echo "Summary: $SUMMARY_FILE" >> "$MASTER_FILE"

echo ""
echo "=== Consolidation Complete ==="
echo "Results saved to: $CONSOLIDATION_DIR"
echo "Summary file: $SUMMARY_FILE"
echo "Master file: $MASTER_FILE"
echo "Pair analyses: $PAIRS_DIR/"
echo ""
echo "Next steps:"
echo "1. Review the pair analyses in $PAIRS_DIR/"
echo "2. Use TOMTOM to compare motifs between pairs"
echo "3. Run the Python consolidation scripts from the pipeline"
