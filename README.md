# STREME Analysis Pipeline# STREME Analysis Pipeline# STREME Analysis Pipeline



A comprehensive, production-ready toolkit for analyzing STREME motif discovery results and their effects on gene expression across multiple genetic lines of *Mimulus guttatus*.



## 🎯 OverviewA comprehensive, production-ready toolkit for analyzing STREME motif discovery results and their effects on gene expression across multiple genetic lines of *Mimulus guttatus*.A comprehensive, production-ready toolkit for analyzing STREME motif discovery results and their effects on gene expression across multiple genetic lines of *Mimulus guttatus*.



This pipeline consolidates STREME motif outputs across multiple genetic lines, validates motif clustering, and provides **mixed-effects statistical analysis** (the gold standard approach in regulatory genomics) to understand motif-expression relationships.



### ✨ Key Features[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)## 🎯 Purpose



- 🔬 **STREME motif consolidation** across genetic lines with IUPAC-aware similarity[![License](https://img.shields.io/badge/license-MIT-green.svg)](LICENSE)

- 📊 **Mixed-effects modeling** - the gold standard in regulatory genomics

- 🎛️ **Multiple statistical approaches** - compare gene-by-gene vs population-level effectsThis pipeline consolidates STREME motif outputs across multiple genetic lines, validates motif clustering, and provides **mixed-effects statistical analysis** (the gold standard approach in regulatory genomics) to understand motif-expression relationships.

- 📈 **Comprehensive analysis** - presence, position bias, and motif variation effects

- 🚀 **Production-ready pipeline** - unified CLI with modular components## 🎯 Overview

- 📚 **Extensive documentation** - guides for setup, usage, and interpretation

## ✨ Key Features

## 🚀 Quick Start

This pipeline consolidates STREME motif outputs across multiple genetic lines, validates motif clustering, and provides **mixed-effects statistical analysis** (the gold standard approach in regulatory genomics) to understand motif-expression relationships.

```bash

# Complete analysis pipeline (recommended)- 🔬 **STREME motif consolidation** across genetic lines with IUPAC-aware similarity

./bin/streme-parser full --input-dir streme_results/ --fasta sequences.fa --expression expr.tsv --output results/

### ✨ Key Features- 📊 **Mixed-effects modeling** - the gold standard in regulatory genomics

# Mixed-effects analysis (if you have features already)

./bin/streme-parser analyze-mixed motif_features.tsv expression.tsv --output analysis/- 🎛️ **Multiple statistical approaches** - compare gene-by-gene vs population-level effects



# Get help on any command- 🔬 **STREME motif consolidation** across genetic lines with IUPAC-aware similarity- 📈 **Comprehensive analysis** - presence, position bias, and motif variation effects

./bin/streme-parser --help

./bin/streme-parser consolidate --help- 📊 **Mixed-effects modeling** - the gold standard in regulatory genomics- 🚀 **Production-ready pipeline** - unified CLI with modular components

```

- 🎛️ **Multiple statistical approaches** - compare gene-by-gene vs population-level effects- 📚 **Extensive documentation** - guides for setup, usage, and interpretation

## 📁 Project Structure

- 📈 **Comprehensive analysis** - presence, position bias, and motif variation effects

```

STREME_parser/- 🚀 **Production-ready pipeline** - unified CLI with modular components## 🚀 Quick Start

├── bin/                    # Main executables (streme-parser CLI)

├── cli_tools/              # Core analysis tools- 📚 **Extensive documentation** - guides for setup, usage, and interpretation

├── pipelines/              # Main orchestration (streme_pipeline.py)

├── docs/                   # Comprehensive documentation```bash

├── scripts/                # Utility scripts

├── archive/                # Alternative implementations## 🚀 Quick Start# Complete analysis pipeline (recommended)

└── requirements.txt        # Python dependencies

```./bin/streme-parser full --input-dir streme_results/ --fasta sequences.fa --expression expr.tsv --output results/



📖 **See [`docs/PROJECT_STRUCTURE.md`](docs/PROJECT_STRUCTURE.md) for detailed structure overview.**```bash



## 📊 Available Commands# Complete analysis pipeline (recommended)# Mixed-effects analysis (if you have features already)



### Core Analysis./bin/streme-parser full \./bin/streme-parser analyze-mixed motif_features.tsv expression.tsv --output analysis/

- **`consolidate`** - Consolidate STREME motifs across genetic lines

- **`validate`** - Validate motif consolidation quality      --input-dir streme_results/ \

- **`extract-features`** - Extract ML features from consolidated motifs

- **`analyze-mixed`** - Mixed-effects analysis (RECOMMENDED approach)    --fasta sequences.fa \# Get help on any command

- **`analyze-comprehensive`** - Advanced regulatory analysis (presence, position, variation)

- **`full`** - Complete pipeline workflow    --expression expr.tsv \./bin/streme-parser --help



### Statistical Approaches    --output results/./bin/streme-parser <command> --help



The pipeline implements multiple modeling approaches:```



1. **Mixed-Effects Models** (recommended)# Mixed-effects analysis (if you have features already)

   - Handles gene-specific baseline differences (random intercepts)

   - Allows motif effects to vary by gene (random slopes)  ./bin/streme-parser analyze-mixed motif_features.tsv expression.tsv --output analysis/## 📁 Project Structure

   - Borrows statistical strength across genes

   - Gold standard in regulatory genomics



2. **Gene-by-Gene Analysis**# Get help on any command```

   - Individual gene modeling

   - No assumptions about shared effects./bin/streme-parser --helpSTREME_parser/

   - Higher interpretability, less statistical power

./bin/streme-parser <command> --help├── bin/                    # Main executables (streme-parser CLI)

3. **Hierarchical Analysis**

   - Groups similar genes for shared effect estimation```├── cli_tools/              # Core analysis tools

   - Balance between power and gene-specificity

├── pipelines/              # Main orchestration (streme_pipeline.py)

## 📚 Documentation

## 📁 Project Structure├── docs/                   # Comprehensive documentation

**Complete guides available in [`docs/`](docs/) directory:**

├── scripts/                # Utility scripts

- **[Getting Started](docs/GETTING_STARTED.md)** - Setup and first steps

- **[Complete Workflow](docs/COMPLETE_WORKFLOW.md)** - Full analysis walkthrough```├── archive/                # Alternative implementations

- **[Statistical Modeling Guide](docs/STATISTICAL_MODELING_GUIDE.md)** - Modeling approaches explained

- **[Project Structure](docs/PROJECT_STRUCTURE.md)** - Detailed project organizationSTREME_parser/└── requirements.txt        # Python dependencies

- **[Output Reference](docs/OUTPUT_COLUMNS_GUIDE.md)** - Output format specifications

├── bin/                    # Main executables (streme-parser CLI)```

## 🎯 Use Cases

├── cli_tools/              # Core analysis tools

- **Regulatory genomics research**: Link motif presence to gene expression differences

- **QTL analysis**: Understand how genetic variants affect regulatory elements├── pipelines/              # Main orchestration (streme_pipeline.py)See [`docs/PROJECT_STRUCTURE.md`](docs/PROJECT_STRUCTURE.md) for detailed structure overview.

- **Comparative genomics**: Analyze motif conservation and divergence across lines

- **Machine learning**: Generate features for predictive modeling of gene expression├── docs/                   # Comprehensive documentation



## 🔧 Installation & Requirements├── scripts/                # Utility scripts## 🚀 Quick Start



```bash├── archive/                # Alternative implementations

# Clone the repository

git clone <repository-url>└── requirements.txt        # Python dependencies### 1. Make the tool executable

cd STREME_parser

``````bash

# Install Python dependencies

pip install -r requirements.txtchmod +x bin/streme-parser



# Make tools executable**📖 See [`docs/PROJECT_STRUCTURE.md`](docs/PROJECT_STRUCTURE.md) for detailed structure overview.**```

chmod +x bin/streme-parser

```



**System Requirements:**## 📊 Available Commands### 2. Consolidate Motifs Across Lines

- Python 3.7+

- pandas, numpy, scikit-learn, matplotlib, seaborn

- statsmodels (for advanced mixed-effects modeling)

### Core Analysis```bash

## 💡 Examples

- **`consolidate`** - Consolidate STREME motifs across genetic lines# Basic consolidation

```bash

# Example 1: Complete pipeline- **`validate`** - Validate motif consolidation quality  ./bin/streme-parser consolidate /path/to/streme/results --output my_analysis

./bin/streme-parser full \

    --input-dir my_streme_results/ \- **`extract-features`** - Extract ML features from consolidated motifs

    --fasta promoter_sequences.fa \

    --expression gene_expression.tsv \- **`analyze-mixed`** - Mixed-effects analysis (RECOMMENDED approach)# With custom parameters

    --output complete_analysis/

- **`analyze-comprehensive`** - Advanced regulatory analysis (presence, position, variation)./bin/streme-parser consolidate /path/to/streme/results \

# Example 2: Mixed-effects analysis only  

./bin/streme-parser analyze-mixed \- **`full`** - Complete pipeline workflow    --output my_analysis/ \

    motif_features.tsv expression_data.tsv \

    --approach all --output mixed_results/    --threshold 0.8 \



# Example 3: Simple feature extraction### Statistical Approaches    --verbose

./bin/streme-parser extract-features \

    consolidated_motifs.tsv sequences.fa \```

    --simple --output simple_features.tsv

```The pipeline implements multiple modeling approaches:



## 🔬 Scientific Background### 3. Validate Consolidation Quality



This pipeline implements statistical approaches commonly used in regulatory genomics:1. **Mixed-Effects Models** (recommended)



- **Mixed-effects models** account for the hierarchical structure of genomic data (observations nested within genes)   - Handles gene-specific baseline differences (random intercepts)```bash

- **IUPAC-aware motif similarity** properly handles degenerate nucleotide codes

- **Multiple testing correction** ensures robust statistical conclusions   - Allows motif effects to vary by gene (random slopes)  # Check clustering quality and overlap handling

- **Cross-validation** provides unbiased performance estimates

   - Borrows statistical strength across genes./bin/streme-parser validate my_analysis/consolidated_streme_sites.tsv

For detailed statistical explanations, see [`docs/STATISTICAL_MODELING_GUIDE.md`](docs/STATISTICAL_MODELING_GUIDE.md).

   - Gold standard in regulatory genomics```

## 🤝 Contributing



This is a research tool developed for *Mimulus guttatus* regulatory analysis. For questions or contributions:

2. **Gene-by-Gene Analysis**### 4. Extract Machine Learning Features

- Open an issue for bug reports or feature requests

- Submit pull requests for improvements   - Individual gene modeling

- Contact the development team for collaboration

   - No assumptions about shared effects```bash

## 📜 License

   - Higher interpretability, less statistical power# Simple binary features (presence/absence only)

[Add license information here]

./bin/streme-parser extract-features my_analysis/consolidated_streme_sites.tsv --simple

## 🙏 Citation

3. **Hierarchical Analysis**

If you use this pipeline in your research, please cite:

   - Groups similar genes for shared effect estimation# Detailed features with expression data

```

[Add citation information here]   - Balance between power and gene-specificity./bin/streme-parser extract-features my_analysis/consolidated_streme_sites.tsv \

```

    --expression expression_data.tsv --top-motifs 100

---

## 📚 Documentation```

🧬 **Transform STREME discoveries into biological insights with statistical rigor** 📊


**Complete guides available in [`docs/`](docs/) directory:**## 📊 Available Commands



- **[Getting Started](docs/GETTING_STARTED.md)** - Setup and first steps### `consolidate` - Consolidate STREME Motifs

- **[Complete Workflow](docs/COMPLETE_WORKFLOW.md)** - Full analysis walkthroughGroups similar motifs across genetic lines using IUPAC-aware similarity scoring.

- **[Statistical Modeling Guide](docs/STATISTICAL_MODELING_GUIDE.md)** - Modeling approaches explained

- **[Project Structure](docs/PROJECT_STRUCTURE.md)** - Detailed project organization### `validate` - Validate Consolidation Quality  

- **[Output Reference](docs/OUTPUT_COLUMNS_GUIDE.md)** - Output format specificationsChecks the quality of motif clustering and identifies potential issues.



## 🎯 Use Cases### `extract-features` - Extract ML Features

Converts consolidated motif data into machine learning-ready feature matrices.

- **Regulatory genomics research**: Link motif presence to gene expression differences- **Simple mode** (`--simple`): Binary presence/absence features only

- **QTL analysis**: Understand how genetic variants affect regulatory elements- **Detailed mode**: Comprehensive features including positional, sequence variation, and density metrics

- **Comparative genomics**: Analyze motif conservation and divergence across lines

- **Machine learning**: Generate features for predictive modeling of gene expression### `full` - Complete Pipeline

Runs consolidation and feature extraction in one command.

## 🔧 Installation & Requirements

## 🛠 Technical Features

```bash

# Clone the repository- **IUPAC-aware motif similarity** with position-by-position comparison

git clone <repository-url>- **Overlap merging** to handle redundant motif hits

cd STREME_parser- **Strand-aware detection** for both forward and reverse complements

- **Position-weighted scoring** (proximity to TSS/gene start matters)

# Install Python dependencies- **Configurable similarity thresholds** for fine-tuning clustering

pip install -r requirements.txt

## � Expected Results

# Make tools executable

chmod +x bin/streme-parser- **Consolidated motif catalog** with unique regulatory elements across all lines

```- **Quality validation reports** showing clustering effectiveness

- **ML-ready feature matrices** for gene expression prediction

**System Requirements:**- **Compatible outputs** for R, Python, and Excel analysis

- Python 3.7+

- pandas, numpy, scikit-learn, matplotlib, seaborn## 🎯 Next Steps

- statsmodels (for advanced mixed-effects modeling)

1. **Run consolidation** on your STREME results

## 💡 Examples2. **Validate** the clustering quality

3. **Extract features** appropriate for your analysis

```bash4. **Correlate** with gene expression data for regulatory insights

# Example 1: Complete pipeline

./bin/streme-parser full \## 📞 Usage Examples

    --input-dir my_streme_results/ \

    --fasta promoter_sequences.fa \```bash

    --expression gene_expression.tsv \# Test on single line

    --output complete_analysis/python cli_tools/gene_motif_mapper.py \

    outputs/consolidated_motif_catalog.txt \

# Example 2: Mixed-effects analysis only      sequences/Genes_IM502_DNA_lifted.fasta \

./bin/streme-parser analyze-mixed \    IM502 \

    motif_features.tsv expression_data.tsv \    --output test_outputs/ \

    --approach all --output mixed_results/    --detailed



# Example 3: Simple feature extraction# Batch process all lines

./bin/streme-parser extract-features \for line in IM502 IM664 IM767 IM1034; do

    consolidated_motifs.tsv sequences.fa \    python cli_tools/gene_motif_mapper.py \

    --simple --output simple_features.tsv        outputs/consolidated_motif_catalog.txt \

```        sequences/Genes_${line}_DNA_lifted.fasta \

        $line \

## 🔬 Scientific Background        --output batch_outputs/

done

This pipeline implements statistical approaches commonly used in regulatory genomics:```



- **Mixed-effects models** account for the hierarchical structure of genomic data (observations nested within genes)---

- **IUPAC-aware motif similarity** properly handles degenerate nucleotide codes

- **Multiple testing correction** ensures robust statistical conclusions**Ready to link regulatory elements to gene expression! 🧬→📊**

- **Cross-validation** provides unbiased performance estimates

For detailed statistical explanations, see [`docs/STATISTICAL_MODELING_GUIDE.md`](docs/STATISTICAL_MODELING_GUIDE.md).

## 🤝 Contributing

This is a research tool developed for *Mimulus guttatus* regulatory analysis. For questions or contributions:

- Open an issue for bug reports or feature requests
- Submit pull requests for improvements
- Contact the development team for collaboration

## 📜 License

[Add license information here]

## 🙏 Citation

If you use this pipeline in your research, please cite:

```
[Add citation information here]
```

---

🧬 **Transform STREME discoveries into biological insights with statistical rigor** 📊