# Project Structure

```
STREME_parser/
├── 📁 bin/                          # Main executables
│   ├── streme-parser                 # Main CLI tool (shell wrapper)
│   └── main.py                       # Python entry point
│
├── 🔧 cli_tools/                     # Core analysis tools
│   ├── motif_consolidator.py         # Consolidate motifs across lines
│   ├── validate_consolidation.py     # Validate motif clustering
│   ├── motif_to_regression_features.py # Extract ML features
│   ├── motif_expression_analyzer.py  # Basic expression analysis
│   ├── mixed_effects_analyzer.py     # Mixed-effects analysis (RECOMMENDED)
│   ├── comprehensive_motif_analyzer.py # Advanced regulatory analysis
│   └── streme_sites_consolidator.py  # Parse STREME sites.tsv files
│
├── 🚀 pipelines/                     # Main orchestration
│   └── streme_pipeline.py            # Master pipeline script
│
├── 📚 docs/                          # All documentation
│   ├── README.md                     # Documentation index
│   ├── GETTING_STARTED.md            # Setup and first steps
│   ├── COMPLETE_WORKFLOW.md          # Full analysis workflow
│   ├── COMPREHENSIVE_REGULATORY_ANALYSIS.md # Advanced analysis guide
│   ├── STATISTICAL_MODELING_GUIDE.md # Statistical approaches explained
│   └── OUTPUT_COLUMNS_GUIDE.md       # Output format reference
│
├── 📜 scripts/                       # Utility scripts
│   ├── consolidate_streme_results.sh # Batch consolidation
│   └── remote_streme_all_lines.sh    # Remote execution helper
│
├── 📦 archive/                       # Alternative implementations
│   ├── advanced_streme_consolidation.py
│   ├── motif_catalog_builder.py
│   └── ...                          # Other experimental tools
│
├── 📋 README.md                      # Main project overview
├── 📋 requirements.txt               # Python dependencies
└── 📋 .gitignore                     # Git ignore patterns
```

## 🎯 Key Components

### Core Tools (`cli_tools/`)
- **motif_consolidator.py**: Groups similar motifs across genetic lines using IUPAC-aware similarity
- **mixed_effects_analyzer.py**: Gold standard statistical analysis for regulatory genomics
- **motif_to_regression_features.py**: Feature extraction for machine learning analysis

### Pipeline (`pipelines/`)
- **streme_pipeline.py**: Unified command-line interface orchestrating all analysis steps

### Documentation (`docs/`)
- Comprehensive guides for setup, usage, and interpretation
- Statistical modeling explanations and best practices
- Output format specifications

### Entry Points (`bin/`)
- **streme-parser**: Main executable providing unified CLI access to all tools
- **main.py**: Python entry point for programmatic access

## 🔄 Analysis Workflow

```mermaid
graph LR
    A[STREME Results] --> B[consolidate]
    B --> C[validate]
    C --> D[extract-features]
    D --> E[analyze-mixed]
    E --> F[Results & Plots]
```

## 🎓 Statistical Approaches

The pipeline implements multiple statistical modeling approaches:

1. **Mixed-Effects Models** (recommended): Gold standard approach
   - Handles gene-specific baseline differences
   - Allows motif effects to vary by gene
   - Borrows statistical strength across genes

2. **Gene-by-Gene Analysis**: Individual gene modeling
   - No assumptions about shared effects
   - Gene-specific effect estimates
   - Higher interpretability but less power

3. **Hierarchical Models**: Groups similar genes
   - Shared effects for similar gene classes
   - Balance between power and specificity

See `docs/STATISTICAL_MODELING_GUIDE.md` for detailed explanations.

## 🚀 Quick Start

```bash
# Complete analysis pipeline
./bin/streme-parser full --input-dir streme_results/ --fasta sequences.fa --expression expr.tsv

# Individual steps
./bin/streme-parser consolidate streme_results/
./bin/streme-parser extract-features motifs.tsv sequences.fa --simple
./bin/streme-parser analyze-mixed features.tsv expression.tsv
```

For detailed usage, see `docs/GETTING_STARTED.md`.