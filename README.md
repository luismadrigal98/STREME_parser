# STREME Analysis Pipeline

A comprehensive toolkit for analyzing STREME motif discovery results across multiple genetic lines of *Mimulus guttatus*.

## 🎯 Purpose

This pipeline consolidates STREME motif outputs across multiple genetic lines, validates motif clustering, and extracts machine learning-ready features for gene expression analysis.

## 📁 Directory Structure

```
STREME_parser/
├── bin/                    # Main executables
│   ├── streme-parser       # Main CLI tool (shell wrapper)
│   └── main.py            # Python entry point
├── cli_tools/              # Core CLI utilities
│   ├── motif_consolidator.py         # Consolidate motifs across lines
│   ├── streme_sites_consolidator.py  # Parse STREME sites.tsv files
│   ├── validate_consolidation.py     # Validate motif clustering quality
│   └── motif_to_regression_features.py # Extract ML features
├── pipelines/              # Main pipeline orchestrator
│   └── streme_pipeline.py  # Master CLI tool
├── scripts/                # Utility scripts for cluster computing
├── outputs/                # Generated results
├── archive/                # Old/deprecated scripts
└── README.md              # This file
```

## 🚀 Quick Start

### 1. Consolidate Motifs Across Lines

```bash
# Basic consolidation
python meme_pipeline.py consolidate /path/to/streme/results

# With custom parameters
python meme_pipeline.py consolidate /path/to/streme/results \
    --output my_outputs/ \
    --threshold 0.8 \
    --verbose
```

**Output:**
- `consolidated_motif_catalog.txt` - Unique motifs with occurrence info
- `regulatory_maps_by_line.txt` - Line-level motif summaries  
- `regulatory_strings_for_ML.fasta` - ML-ready regulatory strings

### 2. Map Motifs to Individual Genes

```bash
# Map motifs to genes for specific line
python meme_pipeline.py map-genes \
    outputs/consolidated_motif_catalog.txt \
    sequences/Genes_IM502_DNA_lifted.fasta \
    IM502 \
    --detailed
```

**Output:**
- `IM502_gene_motif_matrix.txt` - Gene × motif score matrix
- `IM502_detailed_motif_report.txt` - All motif positions per gene

### 3. Full Pipeline (Recommended)

```bash
# Complete analysis for multiple lines
python meme_pipeline.py full \
    /path/to/streme/results \
    /path/to/sequences/ \
    --lines IM502,IM664,IM767,IM1034 \
    --output final_results/
```

## 📊 Key Improvements Over Previous Versions

### ✅ Proper Motif Clustering
- **Before:** `AAAAAAAATTA` and `CTYTCTCTCTCTH` incorrectly clustered together
- **After:** Biologically meaningful clustering with IUPAC-aware similarity

### ✅ Gene-Level Resolution  
- **Before:** Line-level regulatory maps (too coarse)
- **After:** Individual gene motif scores and positions

### ✅ Expression-Ready Outputs
- Gene × motif matrices ready for correlation with expression data
- Position-weighted scoring (TSS proximity matters)
- Strand-aware motif detection

## 🔬 Biological Workflow

```
STREME Results → Motif Consolidation → Gene Mapping → Expression Analysis
     ↓                    ↓                ↓              ↓
Multiple lines    Unique regulatory   Gene-specific   Link regulatory
each with 50-200     elements         motif scores    elements to
motifs per line    across all lines   and positions   gene expression
```

## 📈 Expected Results

- **~1,100 unique motifs** across 10 lines (realistic diversity)
- **~355 common motifs** (shared regulatory elements)  
- **~746 line-specific motifs** (lineage-specific regulation)

## 🛠 Technical Details

### Motif Similarity Algorithm
- IUPAC code expansion and position-by-position comparison
- Length difference penalties (motifs >3bp different won't cluster)
- Sliding window alignment for optimal matching
- Default threshold: 0.75 (adjustable)

### Gene Mapping Features
- Regex-based motif searching with IUPAC support
- Reverse complement detection
- Position-weighted scoring (closer to 5' = higher score)
- Multiple match reporting per gene

### Output Formats
- **TSV matrices** - Compatible with R, Python, Excel
- **FASTA strings** - Ready for ML pipelines
- **Detailed reports** - All motif positions for manual inspection

## 🔗 Integration with Expression Data

The gene-motif matrices are designed to correlate with expression data:

```python
# Example analysis workflow
import pandas as pd

# Load results
motif_scores = pd.read_csv('IM502_gene_motif_matrix.txt', sep='\t', index_col=0)
expression = pd.read_csv('IM502_expression_data.txt', sep='\t', index_col=0)

# Correlate motif presence with expression
correlations = motif_scores.corrwith(expression, axis=0)
significant_motifs = correlations[abs(correlations) > 0.3]
```

## 🎯 Next Steps

1. **Run the pipeline** on your complete dataset
2. **Correlate** gene-motif matrices with expression data
3. **Identify** regulatory motifs associated with specific expression patterns
4. **Validate** findings with known regulatory elements
5. **Build ML models** using regulatory strings

## 📞 Usage Examples

```bash
# Test on single line
python cli_tools/gene_motif_mapper.py \
    outputs/consolidated_motif_catalog.txt \
    sequences/Genes_IM502_DNA_lifted.fasta \
    IM502 \
    --output test_outputs/ \
    --detailed

# Batch process all lines
for line in IM502 IM664 IM767 IM1034; do
    python cli_tools/gene_motif_mapper.py \
        outputs/consolidated_motif_catalog.txt \
        sequences/Genes_${line}_DNA_lifted.fasta \
        $line \
        --output batch_outputs/
done
```

---

**Ready to link regulatory elements to gene expression! 🧬→📊**
