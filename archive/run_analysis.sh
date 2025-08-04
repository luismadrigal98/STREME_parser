#!/bin/bash
# Quick script to run motif analysis on the cluster
# Usage: ./run_analysis.sh

echo "=== MIMULUS MOTIF ANALYSIS ==="
echo "Running comprehensive motif compilation..."
echo

# Run the analysis with different similarity thresholds
echo "Testing with similarity threshold 0.6 (60%)"
python motif_compiler_comprehensive.py /home/l338m483/scratch/MEME_Test/All_lines 0.6

echo
echo "=== ANALYSIS COMPLETE ==="
echo
echo "Output files generated:"
echo "1. comprehensive_motif_catalog.txt"
echo "2. regulatory_maps_by_line.txt" 
echo "3. regulatory_strings_for_ML.fasta"
echo

echo "Quick summary:"
echo "Lines with motifs:"
grep "^## LINE:" regulatory_maps_by_line.txt | wc -l
echo
echo "Total common motifs found:"
grep -c "^cluster_" comprehensive_motif_catalog.txt | head -1
echo
echo "Regulatory strings for ML:"
grep "^>" regulatory_strings_for_ML.fasta | wc -l
