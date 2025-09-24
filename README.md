# STREME Parser — Motif Consolidation and Expression Analysis

STREME Parser is a small toolkit and pipeline for turning STREME motif discovery outputs into consolidated, analysis-ready data. It clusters similar motifs across lines with IUPAC-aware matching, merges overlapping hits, produces a clean consolidated sites table, derives regression-ready features, and includes tools to relate motifs to gene expression across genetic lines.

Works well for multi-line promoter scans where each line was run separately with STREME and produced `streme_*` folders containing a `sites.tsv`.


## Highlights

- IUPAC-aware motif consolidation across lines with length penalties and overlap merging
- True-consensus computation from observed sequences, not just STREME’s consensus pattern
- Rich consolidated output (`consolidated_streme_sites.tsv`) with per-site fields, relative position info, and cluster metadata
- Feature generation for ML/regression (presence, counts, positions, variation) ready to join with expression data
- Two analysis tracks:
  - Motif→Expression regression with multiple models and visualizations
  - Comprehensive analysis including position bias, within-motif sequence variation, and cross-line variation
- A simple pipeline wrapper (`bin/streme-parser`) to run end-to-end steps


## Requirements

- Python 3.8+
- Recommended: Conda environment
- Python packages (install via `requirements.txt`):
  - biopython, pandas, numpy, scipy, statsmodels, scikit-learn
  - matplotlib, seaborn, pyyaml, joblib, tqdm

Install:

```bash
python -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```


## Repository layout

- `bin/`
  - `streme-parser`: Shell wrapper for the main pipeline
  - `main.py`: Delegates to `pipelines/streme_pipeline.py`
- `cli_tools/`
  - `streme_sites_consolidator.py`: Consolidate/validate/features (primary CLI)
  - `motif_to_regression_features.py`: Convert consolidated sites to features
  - `motif_expression_analyzer.py`: Model motif→expression relationships (ABSOLUTE per-line analysis)
  - `relative_motif_analyzer.py`: Model motif→expression relationships (RELATIVE vs IM767 baseline)
  - `comprehensive_motif_analyzer.py`: Position + sequence variation + cross-line effects
  - `validate_consolidation.py`: Standalone validator (duplicated within consolidator CLI)
- `pipelines/`
  - `streme_pipeline.py`: Pipeline entry (consolidate, validate, full)
- `scripts/`
  - `remote_streme_all_lines.sh`: SLURM array job to run RepeatMasker + STREME per FASTA
  - `consolidate_streme_results.sh`: Example shell consolidation helper
- `docs/`: Additional notes and guides


## Inputs and expected files

STREME results directory containing one subfolder per line, e.g.:

```text
streme_results/
├── streme_IM502/sites.tsv
├── streme_IM664/sites.tsv
└── streme_IM767/sites.tsv
```

Each `sites.tsv` should include (STREME defaults) columns similar to:

- `motif_ID`, `seq_ID`, `site_Start`, `site_End`, `site_Strand`, `site_Score`, `site_Sequence`

Expression files (tab-separated, two formats supported):

**Long format (Gene | Line | Expression):**
```text
Gene       Line   Expression
AT1G01010  IM502  0.37
AT1G01010  IM664  -0.12
...
```

**Wide format (gene | LRTadd | IM62 | IM155 | ... | IM767):**
```text
gene       LRTadd  IM62   IM155   IM502   IM664   IM767
AT1G01010  2.45    0.23   0.15    0.37   -0.12    0.00
...
```

The wide format is required for relative analysis (IM767 serves as baseline with expression = 0).


## Quickstart (most direct path)

1. Consolidate all lines’ STREME sites into a single table:

```bash
python cli_tools/streme_sites_consolidator.py consolidate /path/to/streme_results \
  --output outputs/consolidated_streme_sites
```

This writes:

- `outputs/consolidated_streme_sites.tsv`
- `outputs/consolidated_streme_sites_summary.txt`

2. (Optional) Validate consolidation quality and spot suspicious clusters:

```bash
python cli_tools/streme_sites_consolidator.py validate outputs/consolidated_streme_sites.tsv
```

3. Generate regression-ready features (presence, position, sequence variation):

```bash
python cli_tools/streme_sites_consolidator.py features outputs/consolidated_streme_sites.tsv \
  --expression path/to/expression.tsv \
  --top-motifs 100 \
  --min-sites 10 \
  --output-prefix outputs/motif_regression
```

This writes (by default):

- `outputs/motif_regression_features.csv`
- `outputs/motif_regression_feature_descriptions.txt`
- `outputs/motif_regression_summary.txt`

4. Analyze motif effects on expression (three options):

**Option A: Absolute analysis (per-line analysis)**

Analyzes each line independently. Expression data can be in long format (`Gene | Line | Expression`) or wide format (`gene | LRTadd | IM62 | IM155 | ... | IM767`).

```bash
python cli_tools/motif_expression_analyzer.py \
  outputs/consolidated_streme_sites.tsv \
  path/to/expression.tsv \
  --output outputs/absolute_analysis_results
```

**Option B: Relative analysis (comparing to IM767 baseline)**

Generates features comparing each line to IM767 reference. Requires wide format expression data.

```bash
python cli_tools/relative_motif_analyzer.py \
  outputs/consolidated_streme_sites.tsv \
  path/to/expression_wide.tsv \
  --output outputs/relative_analysis_results
```

**Option C: Comprehensive analysis (position bias + sequence variation + cross-line)**

```bash
python cli_tools/comprehensive_motif_analyzer.py \
  outputs/consolidated_streme_sites.tsv \
    path/to/expression.tsv \
    --output outputs/comprehensive_analysis \
    --layers presence position variation cross_line \
    --interactions --top-motifs 50
  ```

  Outputs include: `feature_importance.tsv`, `analysis_summary.md`.


## Pipeline wrapper (recommended)

The wrapper forwards to `pipelines/streme_pipeline.py` and now includes expression analysis commands.

```bash
# Ensure it's executable
chmod +x bin/streme-parser

# Show help
bin/streme-parser --help

# Individual steps
bin/streme-parser consolidate /path/to/streme_results --output outputs/
bin/streme-parser validate outputs/consolidated_streme_sites.tsv
bin/streme-parser analyze outputs/consolidated_streme_sites.tsv expression.tsv --type relative --output results/
bin/streme-parser analyze outputs/consolidated_streme_sites.tsv expression.tsv --type relative --reference-line IM500 --output results/

# Full pipeline (consolidate → validate → analyze)
bin/streme-parser full /path/to/streme_results expression.tsv --output outputs/ --analysis-type relative
bin/streme-parser full /path/to/streme_results expression.tsv --output outputs/ --analysis-type relative --reference-line IM500
```

**Pipeline commands available:**
- `consolidate`: Consolidate STREME motifs across lines
- `validate`: Validate motif consolidation quality  
- `analyze`: Run motif-expression analysis (absolute or relative)
- `full`: Complete pipeline (all three steps)

**Expression analysis types:**
- `--type absolute`: Per-line analysis (each line analyzed independently)
- `--type relative`: Comparative analysis (each line vs reference baseline) - **recommended**
- `--reference-line`: Set reference line for relative analysis (default: IM767)

**Direct pipeline usage:**
```bash
# Relative analysis (recommended)
python pipelines/streme_pipeline.py analyze consolidated_motifs.tsv expression.tsv --type relative

# Full pipeline with expression analysis
python pipelines/streme_pipeline.py full /path/to/streme_results expression.tsv --analysis-type relative
```

Notes:

- The pipeline now includes integrated expression analysis alongside motif processing
- Both absolute and relative analysis types are supported via the wrapper
- For advanced options (detailed features, top motifs), use the CLI tools directly


## Outputs in detail

From consolidation (`consolidated_streme_sites.tsv`):

- Columns include: `consolidated_motif_id`, `original_motif_id`, `motif_consensus`, `original_streme_consensus`,
  `line`, `gene_id`, `start_pos`, `end_pos`, `strand`, `score`, `sequence`,
  `relative_position`, `total_motifs_in_gene`, `relative_position_fraction`, `cluster_size`, `merged_count`, `length`
- Summary file reports counts per line and top motifs

From features (`*_features.csv`):

- Row = one gene×line combination; columns include per-motif presence, counts, position stats, and sequence variation stats
- Companion `*_feature_descriptions.txt` and `*_summary.txt` describe feature definitions and matrix stats

From motif→expression analysis:

- Cross-validated model performance across linear/regularized/tree models
- Feature importance (per-model and aggregate heatmap)
- Predictions and residuals per gene×line (`predictions.tsv`)

From comprehensive analysis:

- Layer-wise R² contributions (presence, position, variation, cross-line)
- Combined model performance and top features


## HPC workflow (SLURM example)

Run STREME per FASTA via array job (example provided in `scripts/remote_streme_all_lines.sh`). Submit with:

```bash
sbatch scripts/remote_streme_all_lines.sh
```

After STREME completes, consolidate:

```bash
python cli_tools/streme_sites_consolidator.py consolidate /path/to/streme_results \
  --output outputs/consolidated_streme_sites
```


## Tips and troubleshooting

- Use Python 3.8+ and install dependencies from `requirements.txt`.
- The consolidator tolerates missing/NaN fields in STREME `sites.tsv` rows and skips incomplete lines.
- If `pandas` isn’t available, output writing falls back to a manual TSV writer.
- Overlap merging is enabled by default in the consolidator to reduce sliding-window artifacts.
- For the regression analyzer, ensure your feature file delimiter matches its expectation (TSV). If you used the feature generator’s CSV output, convert to TSV (see Quickstart).
- Expression file must have headers: `Gene`, `Line`, `Expression`.


## Citation and license

If you use this in a publication, please cite STREME (MEME Suite) for motif discovery. Project license and citation information for this toolkit are TBD.


## Acknowledgements

Built around MEME Suite STREME outputs. Includes utility scripts for an internal regulatory analysis workflow.
