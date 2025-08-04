#!/bin/bash
# Script to run comprehensive motif analysis on STREME results
# Usage: ./run_motif_analysis.sh

# Set paths
STREME_RESULTS_DIR="/home/l338m483/scratch/MEME_Test/All_lines"
PYTHON_SCRIPTS_DIR="/home/l338m483/scratch/MEME_Test/Python"
OUTPUT_DIR="/home/l338m483/scratch/MEME_Test/Results"

# Create output directory
mkdir -p "$OUTPUT_DIR"
cd "$OUTPUT_DIR"

echo "==================================================="
echo "Comprehensive Motif Analysis for Mimulus Lines"
echo "==================================================="
echo "STREME Results Directory: $STREME_RESULTS_DIR"
echo "Python Scripts Directory: $PYTHON_SCRIPTS_DIR"
echo "Output Directory: $OUTPUT_DIR"
echo ""

# Check if the comprehensive script exists
if [ ! -f "$PYTHON_SCRIPTS_DIR/motif_compiler_comprehensive.py" ]; then
    echo "Error: motif_compiler_comprehensive.py not found in $PYTHON_SCRIPTS_DIR"
    echo "Please make sure you copied all Python scripts to the cluster."
    exit 1
fi

# Run the comprehensive motif analysis
echo "Running comprehensive motif analysis..."
echo "This will:"
echo "1. Parse all STREME outputs from 10 lines"
echo "2. Identify common vs unique motifs across lines"
echo "3. Create regulatory maps for each line"
echo "4. Generate regulatory strings for deep learning"
echo ""

python3 "$PYTHON_SCRIPTS_DIR/motif_compiler_comprehensive.py" "$STREME_RESULTS_DIR" 0.8

# Check if analysis completed successfully
if [ $? -eq 0 ]; then
    echo ""
    echo "==================================================="
    echo "Analysis completed successfully!"
    echo "==================================================="
    echo "Output files generated:"
    echo ""
    ls -la *.txt *.fasta 2>/dev/null
    echo ""
    echo "Files description:"
    echo "1. comprehensive_motif_catalog.txt - All motifs with common vs unique classification"
    echo "2. regulatory_maps_by_line.txt - Positional maps for each line (5' to 3')"
    echo "3. regulatory_strings_for_ML.fasta - Regulatory sequences for deep learning"
    echo ""
    echo "Next steps:"
    echo "- Review the motif catalog to see which motifs are common across lines"
    echo "- Examine regulatory maps to understand motif order differences"
    echo "- Use regulatory strings for machine learning models"
else
    echo "Error: Analysis failed. Check the error messages above."
    exit 1
fi
