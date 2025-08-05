# MEME Analysis Pipeline - Complete Usage Example

## Documentation Files
- **[OUTPUT_COLUMNS_GUIDE.md](OUTPUT_COLUMNS_GUIDE.md)** - Detailed explanation of each column in the consolidated TSV output
- **[MOTIF_ANALYSIS_GUIDE.md](MOTIF_ANALYSIS_GUIDE.md)** - Technical guide for motif analysis methods

## Example Workflow for Mimulus guttatus Analysis

### Prerequisites
Your data should be organized like this:
```
your_project/
├── streme_results/
│   ├── streme_Genes_IM502_DNA_lifted/
│   │   ├── sites.tsv       # <- This is what we need!
│   │   ├── streme.txt
│   │   └── streme.html
│   ├── streme_Genes_IM664_DNA_lifted/
│   │   ├── sites.tsv       # <- Position data per motif
│   │   └── ...
│   └── streme_Genes_IM767_DNA_lifted/
│       ├── sites.tsv       # <- Exact sequences found
│       └── ...
└── gene_sequences/         # Optional for this analysis
    ├── Genes_IM502_DNA_lifted.fasta
    └── ...
```

### Method 1: New STREME Sites Analysis (Recommended!)

#### Use the New Sites-Based Consolidator
```bash
cd /path/to/MEME_related/Python/

# Basic analysis - processes all lines
python3 cli_tools/streme_sites_consolidator.py \
    ../All_lines/ \
    --output comprehensive_motif_analysis \
    --verbose

# Custom similarity threshold
python3 cli_tools/streme_sites_consolidator.py \
    ../All_lines/ \
    --threshold 0.8 \
    --output high_stringency_analysis \
    --verbose

# Process specific lines only
python3 cli_tools/streme_sites_consolidator.py \
    ../All_lines/ \
    --lines IM502,IM664,IM767 \
    --output subset_analysis \
    --verbose

# Merge overlapping motif hits (default: enabled)
python3 cli_tools/streme_sites_consolidator.py \
    ../All_lines/ \
    --merge-overlaps \
    --overlap-threshold 0.5 \
    --output merged_analysis \
    --verbose

# Disable overlap merging (keep all individual hits)
python3 cli_tools/streme_sites_consolidator.py \
    ../All_lines/ \
    --no-merge-overlaps \
    --output unmerged_analysis \
    --verbose
```

**Expected output:**
```
=== STREME Sites Consolidator ===
Input directory: ../All_lines/
Similarity threshold: 0.75
Output prefix: comprehensive_motif_analysis

Found 10 STREME directories to process
  Processing streme_Genes_IM502_DNA_lifted/sites.tsv for line IM502
    Found 287 motif sites
  Processing streme_Genes_IM664_DNA_lifted/sites.tsv for line IM664
    Found 312 motif sites
  [... processing all 10 lines ...]

Total sites collected: 2847
Consolidating motifs using similarity threshold 0.75
Found 1243 unique sequences before consolidation
Consolidated into 156 motif clusters

Results written to:
- comprehensive_motif_analysis.tsv
- comprehensive_motif_analysis_summary.txt
```

### Method 2: Legacy Step-by-Step Analysis

#### Step 1: Consolidate all motifs across lines
```bash
cd /path/to/MEME_related/
python3 main.py consolidate ../your_project/streme_results/ \
    --output analysis_outputs/ \
    --threshold 0.75 \
    --verbose
```

**Expected output:**
```
Found 4 STREME directories
Processing streme_IM502... found 45 motifs
Processing streme_IM664... found 52 motifs
Processing streme_IM767... found 41 motifs
Processing streme_IM1034... found 48 motifs

Clustering similar motifs...
Created 178 unique motif clusters
- 42 common motifs (found in 2+ lines)
- 136 line-specific motifs

Results written to:
- analysis_outputs/consolidated_motif_catalog.txt
- analysis_outputs/regulatory_maps_by_line.txt
- analysis_outputs/regulatory_strings_for_ML.fasta
```

#### Step 2: Map motifs to genes for each line
```bash
# Process IM502
python3 main.py map-genes \
    analysis_outputs/consolidated_motif_catalog.txt \
    ../your_project/gene_sequences/Genes_IM502_DNA_lifted.fasta \
    IM502 \
    --output analysis_outputs/ \
    --detailed

# Process IM664
python3 main.py map-genes \
    analysis_outputs/consolidated_motif_catalog.txt \
    ../your_project/gene_sequences/Genes_IM664_DNA_lifted.fasta \
    IM664 \
    --output analysis_outputs/ \
    --detailed

# Continue for other lines...
```

### Method 2: Full Pipeline (Easier!)

```bash
# Run everything at once
python3 main.py full \
    ../your_project/streme_results/ \
    ../your_project/gene_sequences/ \
    --lines IM502,IM664,IM767,IM1034 \
    --output complete_analysis/
```

### Method 3: Using Individual CLI Tools

If you prefer direct control:

```bash
# Use the CLI tools directly
python3 cli_tools/motif_consolidator.py ../your_project/streme_results/

python3 cli_tools/gene_motif_mapper.py \
    consolidated_motif_catalog.txt \
    ../your_project/gene_sequences/ \
    --line IM502 \
    --mismatch 1 \
    --min-score 0.7
```

## Expected Results

### File Outputs
1. **consolidated_motif_catalog.txt** - Master motif list
2. **regulatory_maps_by_line.txt** - Motif summaries per line  
3. **regulatory_strings_for_ML.fasta** - ML-ready sequences
4. **[LINE]_gene_motif_matrix.txt** - Gene × motif matrices (one per line)
5. **[LINE]_detailed_motif_report.txt** - Detailed position reports

### Typical Results for 4 Lines
- **~180 unique motifs** total across all lines
- **~40-50 common motifs** (shared regulatory elements)
- **~130-140 line-specific motifs** (lineage-specific regulation)
- **Gene matrices** ready for correlation with expression data

## Next Steps: Expression Analysis

```python
# Example: Correlate motifs with expression data
import pandas as pd

# Load gene-motif matrix
motifs = pd.read_csv('IM502_gene_motif_matrix.txt', sep='\t', index_col=0)

# Load your expression data (you provide this)
expression = pd.read_csv('IM502_expression_data.txt', sep='\t', index_col=0)

# Find correlations
correlations = motifs.corrwith(expression, axis=0)
significant = correlations[abs(correlations) > 0.3]

print(f"Found {len(significant)} motifs correlated with expression")
```

## Troubleshooting

### Common Issues:
1. **"No STREME directories found"** - Check your input path
2. **"BioPython not found"** - Run `pip install -r requirements.txt`
3. **"Permission denied"** - Make sure scripts are executable: `chmod +x main.py`

### File Format Requirements:
- **STREME output** must contain `streme.txt` files
- **Gene sequences** must be in FASTA format
- **Line names** should match between STREME dirs and FASTA files

### Performance Notes:
- Consolidation: ~2-5 minutes for 4 lines with 50 motifs each
- Gene mapping: ~5-10 minutes per line with 5,000 genes
- Memory usage: ~500MB-1GB depending on dataset size
