# MEME Analysis Pipeline

A comprehensive toolkit for analyzing STREME motif discovery results across multiple genetic lines of *Mimulus guttatus*.

## Quick Start

The pipeline provides a unified CLI tool that can run different analysis subprograms:

```bash
# Basic usage - consolidate motifs
./bin/streme-parser consolidate /path/to/streme/results/ --output my_analysis

# Extract features for machine learning
./bin/streme-parser extract-features my_analysis/consolidated_streme_sites.tsv --simple

# Or using Python directly
python3 pipelines/streme_pipeline.py consolidate /path/to/streme/results/ --output my_analysis
```

## Installation

1. **Prerequisites**: Python 3.6+ with pandas, numpy
```bash
conda activate mycondaenv  # or your preferred environment
pip install pandas numpy
```

2. **Clone/Download** this repository
3. **Make executable**: `chmod +x bin/streme-parser`

## Available Commands

### `consolidate` - Consolidate STREME Sites

Processes STREME `sites.tsv` files from multiple genetic lines, consolidates similar motifs, and creates comprehensive regulatory maps.

#### Key Features:
- ✅ **Motif Consolidation**: Groups similar motifs using IUPAC-aware sequence comparison
- ✅ **Overlap Merging**: Removes redundant overlapping hits (solves poly-A duplication issues)
- ✅ **Repetitive Sequence Handling**: Proper handling of poly-A, poly-T tracts
- ✅ **Exact Coordinates**: Preserves genomic positions and sequences
- ✅ **Comprehensive Output**: 17 columns with detailed motif information

#### Usage Examples:

```bash
```bash
# Basic consolidation with overlap merging (recommended)
./bin/streme-parser consolidate /path/to/streme/results/ --output comprehensive_analysis

# Custom similarity threshold for stricter motif grouping
./bin/streme-parser consolidate /path/to/streme/results/ 
  --threshold 0.8 
  --output strict_analysis
```

# Process specific genetic lines only
./meme consolidate /path/to/streme/results/ \
  --lines IM502,IM664,IM767 \
  --output subset_analysis

# Custom overlap threshold for merging
./meme consolidate /path/to/streme/results/ \
  --overlap-threshold 0.7 \
  --output custom_overlap

# Disable overlap merging (for debugging)
./meme consolidate /path/to/streme/results/ \
  --no-merge-overlaps \
  --output debug_analysis

# Verbose output for troubleshooting
./meme consolidate /path/to/streme/results/ \
  --verbose \
  --output verbose_analysis
```

#### Input Directory Structure:
```
streme_results/
├── streme_IM502/sites.tsv
├── streme_IM664/sites.tsv
├── streme_IM767/sites.tsv
└── streme_IM1034/sites.tsv
```

#### Output Files:
- `{output_prefix}.tsv` - Consolidated motif sites with 17 columns
- `{output_prefix}_summary.txt` - Statistical summary of results

### `extract-features` - Extract Regression Features

Converts consolidated motif data into machine learning ready features for gene expression prediction.

#### Key Features:
- ✅ **Binary Mode**: Simple presence/absence features for each motif (--simple flag)
- ✅ **Detailed Mode**: Comprehensive features including positional, sequence variation, and density metrics
- ✅ **Expression Integration**: Optional expression data integration
- ✅ **Motif Filtering**: Filter by frequency or select top N motifs
- ✅ **Ready for ML**: Output format compatible with scikit-learn and other ML libraries

#### Usage Examples:

```bash
# Simple binary features (presence/absence only)
./bin/streme-parser extract-features consolidated_streme_sites.tsv --simple

# Detailed features with expression data
./bin/streme-parser extract-features consolidated_streme_sites.tsv \
  --expression expression_data.tsv \
  --top-motifs 100

# Filter motifs by minimum occurrence
./bin/streme-parser extract-features consolidated_streme_sites.tsv \
  --min-sites 50 \
  --output-prefix filtered_features
```

#### Feature Types:
- **Simple Mode (--simple)**: Only binary presence/absence features
- **Detailed Mode**: Positional features, sequence variation, GC content, density metrics, spacing regularity

### Future Commands (Coming Soon)

#### `map-genes` - Gene-Level Regulatory Mapping
```bash
./meme map-genes consolidated_results.tsv --output gene_maps/
```

#### `group-seqs` - Sequence Pattern Grouping
```bash
./meme group-seqs consolidated_results.tsv --output sequence_groups/
```

## Output Columns

The consolidated TSV file contains 17 columns with comprehensive motif information:

| Column | Description |
|--------|-------------|
| `consolidated_motif_id` | Unified motif identifier (MOTIF_001, etc.) |
| `original_motif_id` | Original STREME identifier (1-AAAAA, etc.) |
| `motif_consensus` | Consensus sequence pattern |
| `line` | Genetic line identifier |
| `gene_id` | Gene identifier |
| `start_pos`, `end_pos` | Genomic coordinates |
| `strand` | DNA strand (+ or -) |
| `score` | STREME confidence score |
| `sequence` | Actual sequence found |
| `representative_sequence` | Consensus for this motif cluster |
| `relative_position` | Order within gene (1st, 2nd, etc.) |
| `total_motifs_in_gene` | Total regulatory elements in gene |
| `relative_position_fraction` | Normalized position (0.0-1.0) |
| `cluster_size` | Number of variants consolidated |
| `merged_count` | **NEW**: Number of overlapping hits merged |
| `length` | Motif length in base pairs |

👉 **See [`OUTPUT_COLUMNS_GUIDE.md`](OUTPUT_COLUMNS_GUIDE.md) for detailed explanations**

## Key Improvements

### Overlap Merging (NEW!)
Solves the poly-A duplication problem:

**Before** (overlapping duplicates):
```
MOTIF_001  IM1034  MiIM7v11000019m.g  786  800  AAAATAAAAAAAAAG
MOTIF_001  IM1034  MiIM7v11000019m.g  787  801  AAAAATAAAAAAAAA
MOTIF_001  IM1034  MiIM7v11000019m.g  788  802  AAAAAATAAAAAAAA
```

**After** (merged):
```
MOTIF_001  IM1034  MiIM7v11000019m.g  786  805  AAAATAAAAAAAAAG  merged_count: 5
```

### Improved Similarity Algorithm
- Better handling of repetitive sequences (poly-A, poly-T)
- Reduced length penalties for biological repeats
- Maintains precision for complex motifs

## Documentation

- **[`OUTPUT_COLUMNS_GUIDE.md`](OUTPUT_COLUMNS_GUIDE.md)** - Detailed column descriptions
- **[`USAGE_GUIDE.md`](USAGE_GUIDE.md)** - Complete usage examples and workflows
- **[`MOTIF_ANALYSIS_GUIDE.md`](MOTIF_ANALYSIS_GUIDE.md)** - Technical methodology guide

## Typical Results

For 4 genetic lines with ~1000 genes each:
- **Input**: ~100,000 raw STREME motif sites
- **After consolidation**: ~180 unique motif patterns
- **After overlap merging**: ~80,000 non-redundant sites
- **Processing time**: 2-5 minutes

## Troubleshooting

### Common Issues:

1. **Missing pandas**: `pip install pandas`
2. **Permission denied**: `chmod +x meme meme_pipeline.py`
3. **No sites.tsv files found**: Check directory structure
4. **Empty output**: Use `--verbose` flag for debugging

### Getting Help:

```bash
./meme --help                    # General help
./meme consolidate --help        # Consolidation-specific help
```

## Contributing

The pipeline is designed to be modular and extensible. New analysis tools can be added as subcommands in `meme_pipeline.py`.

## Example Workflow

```bash
# 1. Consolidate all STREME results
./meme consolidate ../All_lines/ \
  --output comprehensive_analysis \
  --verbose

# 2. Examine results
head comprehensive_analysis.tsv
cat comprehensive_analysis_summary.txt

# 3. Find genes with high regulatory complexity
awk -F'\t' '$13 > 10' comprehensive_analysis.tsv | cut -f5 | sort | uniq

# 4. Analyze motif distribution patterns
cut -f1 comprehensive_analysis.tsv | sort | uniq -c | sort -nr | head -20
```

This provides a **clean, unified interface** for all your MEME analysis needs!
