# Project Structure

```
STREME_parser/
├── 📁 bin/                          # Main executables
│   ├── streme-parser                 # Main CLI tool (shell wrapper)
│   └── main.py                       # Python entry point
│
├── 🔧 cli_tools/                     # Core analysis tools
│   ├── genome_prep.py                # Genome prep: extract-promoters, mask, background, run-streme
│   ├── genome_terms.py               # Genome/line terminology helpers (canonical "genome", legacy "line")
│   ├── motif_consolidator.py         # Consolidate motifs across genomes
│   ├── validate_consolidation.py     # Validate motif clustering
│   ├── motif_to_regression_features.py # Extract ML features
│   ├── motif_expression_analyzer.py  # Basic expression analysis
│   ├── mixed_effects_analyzer.py     # Mixed-effects analysis (standalone)
│   ├── comprehensive_motif_analyzer.py # Advanced regulatory analysis
│   └── streme_sites_consolidator.py  # Parse/consolidate STREME sites.tsv files (consolidate/validate/features)
│
├── 🚀 pipelines/                     # Main orchestration
│   ├── streme_pipeline.py            # Master pipeline script
│   └── Pipeline_promotor_discovery.yaml # Genome-generic configuration reference
│
├── 📚 docs/                          # All documentation
│   ├── README.md                     # Documentation index
│   ├── GETTING_STARTED.md            # Setup and first steps
│   ├── COMPLETE_WORKFLOW.md          # Full analysis workflow
│   ├── COMPREHENSIVE_REGULATORY_ANALYSIS.md # Advanced analysis guide
│   ├── METHODS.md                    # Materials and Methods (publication-oriented)
│   └── OUTPUT_COLUMNS_GUIDE.md       # Output format reference
│
├── 📜 scripts/                       # Utility scripts
│   ├── consolidate_streme_results.sh # Batch consolidation
│   └── remote_streme_all_lines.sh    # Remote (SLURM array) prepare helper, one genome per task
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
- **genome_prep.py**: From genome assemblies to STREME results — extract promoters (pure Python), mask, build background, run STREME
- **genome_terms.py**: Shared terminology helpers accepting both canonical "genome" and legacy "line" column names
- **motif_consolidator.py**: Groups similar motifs across genomes using IUPAC-aware similarity
- **mixed_effects_analyzer.py**: Standalone mixed-effects statistical analysis for regulatory genomics
- **motif_to_regression_features.py**: Feature extraction for machine learning analysis

### Pipeline (`pipelines/`)
- **streme_pipeline.py**: Unified command-line interface orchestrating all analysis steps (subcommands: `prepare`, `consolidate`, `validate`, `analyze`, `gene-specific`, `full`)
- **Pipeline_promotor_discovery.yaml**: Genome-generic configuration reference for the promoter-discovery workflow

### Documentation (`docs/`)
- Comprehensive guides for setup, usage, and interpretation
- Publication-oriented Materials and Methods (`METHODS.md`)
- Output format specifications

### Entry Points (`bin/`)
- **streme-parser**: Main executable forwarding to `pipelines/streme_pipeline.py`
- **main.py**: Python entry point for programmatic access

## 🔄 Analysis Workflow

```mermaid
graph LR
    A[Genomes + Annotations] --> P[prepare: promoters → mask → background → STREME]
    P --> B[consolidate]
    B --> C[validate]
    C --> D[features]
    D --> E[analyze]
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

See `docs/METHODS.md` for the detailed statistical methodology.

## 🚀 Quick Start

```bash
# Complete pipeline, preparing genomes first from a manifest
./bin/streme-parser full prepared/ expression.tsv \
    --manifest genomes.tsv --jobs 4 --threads 8 \
    --analysis-type relative --reference-genome IM767

# Individual steps
python3 pipelines/streme_pipeline.py prepare --manifest genomes.tsv --output prepared/
./bin/streme-parser consolidate prepared/
python3 cli_tools/streme_sites_consolidator.py features consolidated_streme_sites.tsv --simple
./bin/streme-parser analyze consolidated_streme_sites.tsv expression.tsv --type relative
```

For detailed usage, see `docs/GETTING_STARTED.md`.