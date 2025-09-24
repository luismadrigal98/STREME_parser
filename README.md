# STREME Parser - Comprehensive Motif Discovery and Analysis Pipeline# STREME Parser - Comprehensive Motif Discovery and Analysis Pipeline# STREME Analysis Pipeline



A comprehensive toolkit for motif discovery and regulatory sequence analysis using STREME (Simple, Thorough, Rapid, Enriched Motif Elicitation) and advanced statistical methods.



## OverviewA comprehensive toolkit for motif discovery and regulatory sequence analysis using STREME (Simple, Thorough, Rapid, Enriched Motif Elicitation) and advanced statistical methods.[![Python 3.7+](https://img.shields.io/badge/python-3.7+-blue.svg)](https://www.python.org/downloads/)



This project provides an end-to-end pipeline for:[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

- Motif discovery using STREME

- Advanced motif consolidation and analysis## Overview

- Regulatory sequence analysis

- Statistical modeling and visualizationA comprehensive, production-ready toolkit for analyzing STREME motif discovery results and their effects on gene expression across multiple genetic lines of *Mimulus guttatus*.

- Mixed-effects modeling for complex experimental designs

This project provides an end-to-end pipeline for:

## Features

- Motif discovery using STREME## 🎯 Overview

### Core Functionality

- **STREME Integration**: Automated motif discovery with customizable parameters- Advanced motif consolidation and analysis

- **Motif Consolidation**: Advanced algorithms to merge similar motifs and reduce redundancy

- **Regulatory Analysis**: Comprehensive analysis of regulatory sequences and their properties- Regulatory sequence analysisThis pipeline consolidates STREME motif outputs across multiple genetic lines, validates motif clustering, and provides **mixed-effects statistical analysis** (the gold standard approach in regulatory genomics) to understand motif-expression relationships.

- **Statistical Modeling**: Mixed-effects models and regression analysis for motif-expression relationships

- **Visualization**: Rich plotting and visualization capabilities for results interpretation- Statistical modeling and visualization



### Key Components- Mixed-effects modeling for complex experimental designs### ✨ Key Features

- **Pipeline Orchestration**: YAML-based configuration for reproducible analyses

- **CLI Tools**: Command-line interfaces for batch processing and automation

- **Modular Architecture**: Extensible design for custom analysis workflows

- **Quality Control**: Built-in validation and consolidation verification tools## Features- 🔬 **STREME motif consolidation** across genetic lines with IUPAC-aware similarity



## Installation- 📊 **Mixed-effects modeling** - the gold standard in regulatory genomics



### Prerequisites### Core Functionality- 🎛️ **Multiple statistical approaches** - compare gene-by-gene vs population-level effects

- Python 3.8 or higher

- MEME Suite (for STREME)- **STREME Integration**: Automated motif discovery with customizable parameters- 📈 **Comprehensive analysis** - presence, position bias, and motif variation effects

- Git

- **Motif Consolidation**: Advanced algorithms to merge similar motifs and reduce redundancy- 🚀 **Production-ready pipeline** - unified CLI with modular components

### Setup

```bash- **Regulatory Analysis**: Comprehensive analysis of regulatory sequences and their properties- 📚 **Extensive documentation** - guides for setup, usage, and interpretation

# Clone the repository

git clone https://github.com/luismadrigal98/STREME_parser.git- **Statistical Modeling**: Mixed-effects models and regression analysis for motif-expression relationships

cd STREME_parser

- **Visualization**: Rich plotting and visualization capabilities for results interpretation## 🚀 Quick Start

# Install dependencies

pip install -r requirements.txt



# Make scripts executable (if needed)### Key Components```bash

chmod +x scripts/*.sh

chmod +x bin/streme-parser- **Pipeline Orchestration**: YAML-based configuration for reproducible analyses# Complete analysis pipeline (recommended)

```

- **CLI Tools**: Command-line interfaces for batch processing and automation./bin/streme-parser full --input-dir streme_results/ --fasta sequences.fa --expression expr.tsv --output results/

## Quick Start

- **Modular Architecture**: Extensible design for custom analysis workflows

### Basic Motif Discovery

```bash- **Quality Control**: Built-in validation and consolidation verification tools# Mixed-effects analysis (if you have features already)

# Run STREME analysis on your sequences

python pipelines/streme_pipeline.py --input sequences.fasta --output results/./bin/streme-parser analyze-mixed motif_features.tsv expression.tsv --output analysis/



# Consolidate and analyze results## Installation

python cli_tools/comprehensive_motif_analyzer.py --input results/ --output analysis/

```# Get help on any command



### Using the Main Interface### Prerequisites./bin/streme-parser --help

```bash

# Interactive mode- Python 3.8 or higher```

python bin/main.py

- MEME Suite (for STREME)

# Command-line mode

python bin/main.py --mode batch --config config.yaml- Git## 📁 Project Structure

```



## Project Structure

### Setup```

```

├── bin/                    # Main executables and entry points```bashSTREME_parser/

├── pipelines/             # Analysis pipelines and workflows  

├── cli_tools/             # Command-line analysis tools# Clone the repository├── bin/                    # Main executables (streme-parser CLI)

├── scripts/               # Utility scripts and automation

├── docs/                  # Comprehensive documentationgit clone https://github.com/luismadrigal98/STREME_parser.git├── cli_tools/              # Core analysis tools

├── outputs/               # Default output directory

├── archive/               # Legacy and experimental codecd STREME_parser├── pipelines/              # Main orchestration (streme_pipeline.py)

└── requirements.txt       # Python dependencies

```├── docs/                   # Comprehensive documentation



## Usage# Install dependencies├── scripts/                # Utility scripts



### Motif Discovery Pipelinepip install -r requirements.txt├── archive/                # Alternative implementations

The STREME pipeline provides automated motif discovery:

└── requirements.txt        # Python dependencies

```python

from pipelines.streme_pipeline import STREMEPipeline# Make scripts executable (if needed)```



pipeline = STREMEPipeline(chmod +x scripts/*.sh

    input_fasta="sequences.fasta",

    output_dir="results/",chmod +x bin/streme-parser📖 **See [`docs/PROJECT_STRUCTURE.md`](docs/PROJECT_STRUCTURE.md) for detailed structure overview.**

    min_width=6,

    max_width=20```

)

pipeline.run()## 📊 Available Commands

```

## Quick Start

### Motif Consolidation

Consolidate similar motifs to reduce redundancy:### Core Analysis



```bash### Basic Motif Discovery

python cli_tools/motif_consolidator.py \

    --input results/streme_output/ \```bash- **`consolidate`** - Consolidate STREME motifs across genetic lines

    --output consolidated/ \

    --similarity-threshold 0.8# Run STREME analysis on your sequences- **`validate`** - Validate motif consolidation quality

```

python pipelines/streme_pipeline.py --input sequences.fasta --output results/- **`extract-features`** - Extract ML features from consolidated motifs

### Regulatory Analysis

Perform comprehensive regulatory sequence analysis:- **`analyze-mixed`** - Mixed-effects analysis (RECOMMENDED approach)



```bash# Consolidate and analyze results- **`analyze-comprehensive`** - Advanced regulatory analysis (presence, position, variation)

python cli_tools/comprehensive_motif_analyzer.py \

    --motifs consolidated/motifs.txt \python cli_tools/comprehensive_motif_analyzer.py --input results/ --output analysis/- **`full`** - Complete pipeline workflow

    --sequences sequences.fasta \

    --expression expression_data.csv \```

    --output analysis/

```### Statistical Approaches



### Mixed-Effects Modeling### Using the Main Interface

For complex experimental designs with multiple factors:

```bash1. **Mixed-Effects Models** (recommended)

```bash

python cli_tools/mixed_effects_analyzer.py \# Interactive mode   - Handles gene-specific baseline differences (random intercepts)

    --data analysis_results.csv \

    --design experimental_design.csv \python bin/main.py   - Allows motif effects to vary by gene (random slopes)

    --output mixed_effects/

```   - Borrows statistical strength across genes



## Configuration# Command-line mode   - Gold standard in regulatory genomics



### Pipeline Configurationpython bin/main.py --mode batch --config config.yaml

Use YAML files for reproducible analysis configurations:

```2. **Gene-by-Gene Analysis**

```yaml

# Pipeline_promotor_discovery.yaml   - Individual gene modeling

streme:

  min_width: 6## Project Structure   - No assumptions about shared effects

  max_width: 20

  evt: 0.05   - Higher interpretability, less statistical power

  

consolidation:```

  similarity_threshold: 0.8

  min_sites: 5├── bin/                    # Main executables and entry points3. **Hierarchical Analysis**

  

analysis:├── pipelines/             # Analysis pipelines and workflows     - Groups similar genes for shared effect estimation

  statistical_tests: true

  generate_plots: true├── cli_tools/             # Command-line analysis tools   - Balance between power and gene-specificity

```

├── scripts/               # Utility scripts and automation

### Environment Variables

Set these environment variables for optimal performance:├── docs/                  # Comprehensive documentation## 🔧 Installation & Requirements



```bash├── outputs/               # Default output directory

export MEME_BIN_PATH=/path/to/meme/bin

export STREME_PARSER_THREADS=4├── archive/               # Legacy and experimental code```bash

export STREME_PARSER_MEMORY=8G

```└── requirements.txt       # Python dependencies# Clone the repository



## Output Files```git clone <repository-url>



### Key Output Filescd STREME_parser

- `motifs_consolidated.txt` - Consolidated motif definitions

- `site_analysis.csv` - Detailed binding site analysis  ## Usage

- `regulatory_features.csv` - Extracted regulatory features

- `expression_correlations.csv` - Motif-expression relationships# Install Python dependencies

- `statistical_summary.txt` - Analysis summary and statistics

### 1. Motif Discovery Pipelinepip install -r requirements.txt

### Visualization Outputs

- Motif logos and sequence logosThe STREME pipeline (`pipelines/streme_pipeline.py`) provides automated motif discovery:

- Correlation heatmaps

- Distribution plots# Make tools executable

- Statistical model diagnostics

```pythonchmod +x bin/streme-parser

## Advanced Usage

from pipelines.streme_pipeline import STREMEPipeline```

### Custom Analysis Workflows

Extend the pipeline with custom analysis modules:



```pythonpipeline = STREMEPipeline(**System Requirements:**

from pipelines.streme_pipeline import STREMEPipeline

from cli_tools.comprehensive_motif_analyzer import MotifAnalyzer    input_fasta="sequences.fasta",- Python 3.7+



# Custom workflow    output_dir="results/",- pandas, numpy, scikit-learn, matplotlib, seaborn, statsmodels

pipeline = STREMEPipeline(config='custom_config.yaml')

results = pipeline.run()    min_width=6,



analyzer = MotifAnalyzer()    max_width=20

analyzer.load_results(results)

analyzer.add_custom_analysis(my_custom_function))

final_results = analyzer.run()

```pipeline.run()## 💡 Examples## 📊 Available Commands# Complete analysis pipeline (recommended)# Mixed-effects analysis (if you have features already)



### Batch Processing```

Process multiple datasets efficiently:



```bash

# Batch analysis script### 2. Motif Consolidation

./scripts/remote_streme_all_lines.sh input_directory/ output_directory/

```Consolidate similar motifs to reduce redundancy:```bash



## Documentation



Comprehensive documentation is available in the `docs/` directory:```bash# Example 1: Complete pipeline



- [`GETTING_STARTED.md`](docs/GETTING_STARTED.md) - Detailed setup and first stepspython cli_tools/motif_consolidator.py \

- [`PROJECT_STRUCTURE.md`](docs/PROJECT_STRUCTURE.md) - Architecture and organization

- [`COMPLETE_WORKFLOW.md`](docs/COMPLETE_WORKFLOW.md) - End-to-end analysis guide    --input results/streme_output/ \./bin/streme-parser full \### Core Analysis./bin/streme-parser full \./bin/streme-parser analyze-mixed motif_features.tsv expression.tsv --output analysis/

- [`COMPREHENSIVE_REGULATORY_ANALYSIS.md`](docs/COMPREHENSIVE_REGULATORY_ANALYSIS.md) - Advanced analysis methods

- [`OUTPUT_COLUMNS_GUIDE.md`](docs/OUTPUT_COLUMNS_GUIDE.md) - Output file specifications    --output consolidated/ \



## Requirements    --similarity-threshold 0.8    --input-dir my_streme_results/ \



- Python 3.8+```

- NumPy >= 1.19.0

- Pandas >= 1.3.0    --fasta promoter_sequences.fa \- **`consolidate`** - Consolidate STREME motifs across genetic lines

- Biopython >= 1.78

- SciPy >= 1.7.0### 3. Regulatory Analysis

- Scikit-learn >= 1.0.0

- Statsmodels >= 0.12.0Perform comprehensive regulatory sequence analysis:    --expression gene_expression.tsv \

- Matplotlib >= 3.3.0

- Seaborn >= 0.11.0

- PyYAML >= 5.4.0

```bash    --output complete_analysis/- **`validate`** - Validate motif consolidation quality      --input-dir streme_results/ \

See [`requirements.txt`](requirements.txt) for complete dependency list.

python cli_tools/comprehensive_motif_analyzer.py \

## Contributing

    --motifs consolidated/motifs.txt \

1. Fork the repository

2. Create a feature branch (`git checkout -b feature/new-analysis`)    --sequences sequences.fasta \

3. Make your changes with appropriate tests

4. Submit a pull request with detailed description    --expression expression_data.csv \# Example 2: Mixed-effects analysis only- **`extract-features`** - Extract ML features from consolidated motifs



## License    --output analysis/



This project is licensed under the MIT License - see the LICENSE file for details.```./bin/streme-parser analyze-mixed \



## Citation



If you use this toolkit in your research, please cite:### 4. Mixed-Effects Modeling    motif_features.tsv expression_data.tsv \- **`analyze-mixed`** - Mixed-effects analysis (RECOMMENDED approach)    --fasta sequences.fa \# Get help on any command



```For complex experimental designs with multiple factors:

STREME Parser: Comprehensive Motif Discovery and Regulatory Analysis Pipeline

[Your Name et al.]    --approach all --output mixed_results/

[Year]

``````bash



## Supportpython cli_tools/mixed_effects_analyzer.py \- **`analyze-comprehensive`** - Advanced regulatory analysis (presence, position, variation)



- **Issues**: Report bugs and feature requests via GitHub Issues    --data analysis_results.csv \

- **Documentation**: Check the `docs/` directory for detailed guides

- **Contact**: [Your contact information]    --design experimental_design.csv \# Example 3: Simple feature extraction



## Changelog    --output mixed_effects/



### Version 2.0.0 (Current)```./bin/streme-parser extract-features \- **`full`** - Complete pipeline workflow    --expression expr.tsv \./bin/streme-parser --help

- Added comprehensive motif consolidation algorithms

- Implemented mixed-effects modeling capabilities

- Enhanced visualization and reporting features

- Improved pipeline configurability and modularity## Configuration    consolidated_motifs.tsv sequences.fa \



### Version 1.0.0

- Initial release with basic STREME integration

- Core motif discovery and analysis functionality### Pipeline Configuration    --simple --output simple_features.tsv

Use YAML files for reproducible analysis configurations:

```

```yaml

# Pipeline_promotor_discovery.yaml### Statistical Approaches    --output results/./bin/streme-parser <command> --help

streme:

  min_width: 6## 📚 Documentation

  max_width: 20

  evt: 0.05

  

consolidation:**Complete guides available in [`docs/`](docs/) directory:**

  similarity_threshold: 0.8

  min_sites: 5The pipeline implements multiple modeling approaches:```

  

analysis:- **[Getting Started](docs/GETTING_STARTED.md)** - Setup and first steps

  statistical_tests: true

  generate_plots: true- **[Complete Workflow](docs/COMPLETE_WORKFLOW.md)** - Full analysis walkthrough

```

- **[Statistical Modeling Guide](docs/STATISTICAL_MODELING_GUIDE.md)** - Modeling approaches explained

### Environment Variables

Set these environment variables for optimal performance:- **[Project Structure](docs/PROJECT_STRUCTURE.md)** - Detailed project organization1. **Mixed-Effects Models** (recommended)# Mixed-effects analysis (if you have features already)



```bash- **[Output Reference](docs/OUTPUT_COLUMNS_GUIDE.md)** - Output format specifications

export MEME_BIN_PATH=/path/to/meme/bin

export STREME_PARSER_THREADS=4   - Handles gene-specific baseline differences (random intercepts)

export STREME_PARSER_MEMORY=8G

```## 🎯 Use Cases



## Output Files   - Allows motif effects to vary by gene (random slopes)  ./bin/streme-parser analyze-mixed motif_features.tsv expression.tsv --output analysis/## 📁 Project Structure



### Key Output Files- **Regulatory genomics research**: Link motif presence to gene expression differences

- `motifs_consolidated.txt` - Consolidated motif definitions

- `site_analysis.csv` - Detailed binding site analysis  - **QTL analysis**: Understand how genetic variants affect regulatory elements   - Borrows statistical strength across genes

- `regulatory_features.csv` - Extracted regulatory features

- `expression_correlations.csv` - Motif-expression relationships- **Comparative genomics**: Analyze motif conservation and divergence across lines

- `statistical_summary.txt` - Analysis summary and statistics

- **Machine learning**: Generate features for predictive modeling of gene expression   - Gold standard in regulatory genomics

### Visualization Outputs

- Motif logos and sequence logos

- Correlation heatmaps

- Distribution plots## 🔬 Scientific Background

- Statistical model diagnostics



## Advanced Usage

This pipeline implements statistical approaches commonly used in regulatory genomics:2. **Gene-by-Gene Analysis**# Get help on any command```

### Custom Analysis Workflows

Extend the pipeline with custom analysis modules:



```python- **Mixed-effects models** account for the hierarchical structure of genomic data (observations nested within genes)   - Individual gene modeling

from pipelines.streme_pipeline import STREMEPipeline

from cli_tools.comprehensive_motif_analyzer import MotifAnalyzer- **IUPAC-aware motif similarity** properly handles degenerate nucleotide codes



# Custom workflow- **Multiple testing correction** ensures robust statistical conclusions   - No assumptions about shared effects./bin/streme-parser --helpSTREME_parser/

pipeline = STREMEPipeline(config='custom_config.yaml')

results = pipeline.run()- **Cross-validation** provides unbiased performance estimates



analyzer = MotifAnalyzer()   - Higher interpretability, less statistical power

analyzer.load_results(results)

analyzer.add_custom_analysis(my_custom_function)For detailed statistical explanations, see [`docs/STATISTICAL_MODELING_GUIDE.md`](docs/STATISTICAL_MODELING_GUIDE.md).

final_results = analyzer.run()

```./bin/streme-parser <command> --help├── bin/                    # Main executables (streme-parser CLI)



### Batch Processing## 🤝 Contributing

Process multiple datasets efficiently:

3. **Hierarchical Analysis**

```bash

# Batch analysis scriptThis is a research tool developed for *Mimulus guttatus* regulatory analysis. For questions or contributions:

./scripts/remote_streme_all_lines.sh input_directory/ output_directory/

```   - Groups similar genes for shared effect estimation```├── cli_tools/              # Core analysis tools



## Documentation- Open an issue for bug reports or feature requests



Comprehensive documentation is available in the `docs/` directory:- Submit pull requests for improvements   - Balance between power and gene-specificity



- [`GETTING_STARTED.md`](docs/GETTING_STARTED.md) - Detailed setup and first steps- Contact the development team for collaboration

- [`PROJECT_STRUCTURE.md`](docs/PROJECT_STRUCTURE.md) - Architecture and organization

- [`COMPLETE_WORKFLOW.md`](docs/COMPLETE_WORKFLOW.md) - End-to-end analysis guide├── pipelines/              # Main orchestration (streme_pipeline.py)

- [`COMPREHENSIVE_REGULATORY_ANALYSIS.md`](docs/COMPREHENSIVE_REGULATORY_ANALYSIS.md) - Advanced analysis methods

- [`OUTPUT_COLUMNS_GUIDE.md`](docs/OUTPUT_COLUMNS_GUIDE.md) - Output file specifications## 📜 License



## Requirements## 📚 Documentation



- Python 3.8+This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

- NumPy >= 1.19.0

- Pandas >= 1.3.0## 📁 Project Structure├── docs/                   # Comprehensive documentation

- Biopython >= 1.78

- SciPy >= 1.7.0## 🙏 Citation

- Scikit-learn >= 1.0.0

- Statsmodels >= 0.12.0**Complete guides available in [`docs/`](docs/) directory:**

- Matplotlib >= 3.3.0

- Seaborn >= 0.11.0If you use this pipeline in your research, please cite:

- PyYAML >= 5.4.0

├── scripts/                # Utility scripts

See [`requirements.txt`](requirements.txt) for complete dependency list.

```

## Contributing

[Add citation information here]- **[Getting Started](docs/GETTING_STARTED.md)** - Setup and first steps

1. Fork the repository

2. Create a feature branch (`git checkout -b feature/new-analysis`)```

3. Make your changes with appropriate tests

4. Submit a pull request with detailed description- **[Complete Workflow](docs/COMPLETE_WORKFLOW.md)** - Full analysis walkthrough```├── archive/                # Alternative implementations



## License---



This project is licensed under the MIT License - see the LICENSE file for details.- **[Statistical Modeling Guide](docs/STATISTICAL_MODELING_GUIDE.md)** - Modeling approaches explained



## Citation🧬 **Transform STREME discoveries into biological insights with statistical rigor** 📊

- **[Project Structure](docs/PROJECT_STRUCTURE.md)** - Detailed project organizationSTREME_parser/└── requirements.txt        # Python dependencies

If you use this toolkit in your research, please cite:

- **[Output Reference](docs/OUTPUT_COLUMNS_GUIDE.md)** - Output format specifications

```

STREME Parser: Comprehensive Motif Discovery and Regulatory Analysis Pipeline├── bin/                    # Main executables (streme-parser CLI)```

[Your Name et al.]

[Year]## 🎯 Use Cases

```

├── cli_tools/              # Core analysis tools

## Support

- **Regulatory genomics research**: Link motif presence to gene expression differences

- **Issues**: Report bugs and feature requests via GitHub Issues

- **Documentation**: Check the `docs/` directory for detailed guides- **QTL analysis**: Understand how genetic variants affect regulatory elements├── pipelines/              # Main orchestration (streme_pipeline.py)See [`docs/PROJECT_STRUCTURE.md`](docs/PROJECT_STRUCTURE.md) for detailed structure overview.

- **Contact**: [Your contact information]

- **Comparative genomics**: Analyze motif conservation and divergence across lines

## Changelog

- **Machine learning**: Generate features for predictive modeling of gene expression├── docs/                   # Comprehensive documentation

### Version 2.0.0 (Current)

- Added comprehensive motif consolidation algorithms

- Implemented mixed-effects modeling capabilities

- Enhanced visualization and reporting features## 🔧 Installation & Requirements├── scripts/                # Utility scripts## 🚀 Quick Start

- Improved pipeline configurability and modularity



### Version 1.0.0

- Initial release with basic STREME integration```bash├── archive/                # Alternative implementations

- Core motif discovery and analysis functionality
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