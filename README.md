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

### 1. Make the tool executable
```bash
chmod +x bin/streme-parser
```

### 2. Consolidate Motifs Across Lines

```bash
# Basic consolidation
./bin/streme-parser consolidate /path/to/streme/results --output my_analysis

# With custom parameters
./bin/streme-parser consolidate /path/to/streme/results \
    --output my_analysis/ \
    --threshold 0.8 \
    --verbose
```

### 3. Validate Consolidation Quality

```bash
# Check clustering quality and overlap handling
./bin/streme-parser validate my_analysis/consolidated_streme_sites.tsv
```

### 4. Extract Machine Learning Features

```bash
# Simple binary features (presence/absence only)
./bin/streme-parser extract-features my_analysis/consolidated_streme_sites.tsv --simple

# Detailed features with expression data
./bin/streme-parser extract-features my_analysis/consolidated_streme_sites.tsv \
    --expression expression_data.tsv --top-motifs 100
```

## 📊 Available Commands

### `consolidate` - Consolidate STREME Motifs
Groups similar motifs across genetic lines using IUPAC-aware similarity scoring.

### `validate` - Validate Consolidation Quality  
Checks the quality of motif clustering and identifies potential issues.

### `extract-features` - Extract ML Features
Converts consolidated motif data into machine learning-ready feature matrices.
- **Simple mode** (`--simple`): Binary presence/absence features only
- **Detailed mode**: Comprehensive features including positional, sequence variation, and density metrics

### `full` - Complete Pipeline
Runs consolidation and feature extraction in one command.

## 🛠 Technical Features

- **IUPAC-aware motif similarity** with position-by-position comparison
- **Overlap merging** to handle redundant motif hits
- **Strand-aware detection** for both forward and reverse complements
- **Position-weighted scoring** (proximity to TSS/gene start matters)
- **Configurable similarity thresholds** for fine-tuning clustering

## � Expected Results

- **Consolidated motif catalog** with unique regulatory elements across all lines
- **Quality validation reports** showing clustering effectiveness
- **ML-ready feature matrices** for gene expression prediction
- **Compatible outputs** for R, Python, and Excel analysis

## 🎯 Next Steps

1. **Run consolidation** on your STREME results
2. **Validate** the clustering quality
3. **Extract features** appropriate for your analysis
4. **Correlate** with gene expression data for regulatory insights

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
