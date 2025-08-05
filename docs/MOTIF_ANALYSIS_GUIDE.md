# Motif Analysis Pipeline - Usage Guide

## Overview
This pipeline helps you compile and analyze motifs from independent STREME runs across multiple lines, identifying common vs unique motifs for understanding expression differences and deep learning applications.

## Files Created:

1. **motif_compiler_comprehensive.py** - Main analysis script
2. **consolidate_streme_output.py** - Simple motif extraction script
3. **group_sequences_simple.py** - Basic sequence grouping script

## Usage

### Step 1: Compile All Motifs from Multiple Lines

```bash
python motif_compiler_comprehensive.py /path/to/streme_results [similarity_threshold]
```

**Arguments:**
- `streme_results_directory`: Directory containing your STREME output folders (e.g., streme_line1, streme_line2, etc.)
- `similarity_threshold`: Optional, default 0.8 (80% similarity to group motifs)

**Example:**
```bash
python motif_compiler_comprehensive.py /home/l338m483/scratch/MEME_Test/All_lines 0.8
```

This will generate:
- `comprehensive_motif_catalog.txt` - All motifs with common vs unique classification
- `regulatory_maps_by_line.txt` - Positional maps for each line
- `regulatory_strings_for_ML.fasta` - Regulatory sequences for deep learning

### Step 2: Analyze Results

#### comprehensive_motif_catalog.txt
Shows:
- Summary statistics of total clusters, common motifs, line-specific motifs
- List of common motifs (present in multiple lines)
- List of line-specific motifs
- Detailed cluster information with consensus sequences

#### regulatory_maps_by_line.txt
For each line shows:
- Total motifs found
- Number of common vs unique motifs
- Regulatory string (motif order from 5' to 3')
- Detailed position table

#### regulatory_strings_for_ML.fasta
FASTA format sequences where:
- Each line has a regulatory "sequence" 
- Motif clusters are represented as elements
- Can be used directly in ML/DL models

## Key Features

### Motif Comparison Across Lines
- **Common motifs**: Found in multiple lines (may indicate conserved regulatory elements)
- **Unique motifs**: Found only in specific lines (may explain expression differences)
- **Similarity clustering**: Groups similar motifs using consensus sequence comparison

### Regulatory String Generation
- Creates compressed representations of regulatory landscapes
- Maintains 5' to 3' order information
- Suitable for machine learning input

### Deep Learning Applications
- Regulatory strings can be treated as sequences for:
  - RNN/LSTM models
  - Transformer models
  - Convolutional neural networks
  - Expression prediction models

## Understanding Your Results

### Expression Differences Analysis
1. **Common motifs** suggest shared regulatory mechanisms
2. **Line-specific motifs** may explain expression differences between lines
3. **Motif order** (regulatory strings) captures combinatorial effects

### Machine Learning Pipeline
```python
# Example: Load regulatory strings for ML
sequences = []
labels = []

with open('regulatory_strings_for_ML.fasta', 'r') as f:
    for line in f:
        if line.startswith('>'):
            label = line.strip()[1:]  # Line name
            labels.append(label)
        else:
            sequences.append(line.strip())

# Now you can encode sequences for your ML model
```

## Advanced Options

### Adjusting Similarity Threshold
- **Higher threshold (0.9-1.0)**: More stringent, fewer common motifs
- **Lower threshold (0.6-0.8)**: More permissive, more common motifs
- **Default (0.8)**: Good balance for most applications

### Custom Analysis
You can modify the script to:
- Use PWM-based similarity instead of consensus sequences
- Add expression data integration
- Customize regulatory string encoding
- Add motif functional annotation

## Next Steps for Deep Learning

1. **Encode regulatory strings** as numerical vectors
2. **Add expression data** as target variables
3. **Train models** to predict expression from regulatory patterns
4. **Identify important motif combinations** using model interpretation

This approach will help you understand how motif presence and order impact gene expression differences between your Mimulus lines.
