# STREME Parser — Promoter Discovery, Motif Consolidation, and Expression Analysis

STREME Parser is a toolkit and pipeline for regulatory-element discovery across **genomes**. Starting from genome assemblies plus annotations, it extracts candidate promoters (upstream-of-TSS regions), masks them, runs STREME motif discovery, then clusters similar motifs across genomes with IUPAC-aware matching, merges overlapping hits, produces a clean consolidated sites table, derives regression-ready features, and relates motifs to gene expression.

It is genome-generic: use it for several genomes of a genus (e.g. *Penstemon virgatus*, *P. barbatus*) or for multiple ecotypes/inbred lines of one species (e.g. *Mimulus guttatus* IM lines). Provide one or many `(genome, annotation)` pairs — multiple genomes are processed in parallel; a single genome parallelises its heavy steps.

> **Terminology:** the canonical term is **genome**. Earlier versions used "line"; the legacy `line`/`Line` column names and the `--reference-line` flag are still accepted everywhere, so existing Mimulus data needs no migration.


## Highlights

- IUPAC-aware motif consolidation across genomes with length penalties and overlap merging
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
  - `genome_prep.py`: Genome preparation — `extract-promoters` (pure Python), `mask`, `background`, `run-streme`
  - `genome_terms.py`: Shared genome/line terminology helpers (canonical "genome", legacy "line" alias)
  - `streme_sites_consolidator.py`: Consolidate/validate/features (primary CLI)
  - `motif_to_regression_features.py`: Convert consolidated sites to features
  - `motif_expression_analyzer.py`: Model motif→expression relationships (ABSOLUTE per-genome analysis)
  - `relative_motif_analyzer.py`: Model motif→expression relationships (RELATIVE vs a reference genome)
  - `comprehensive_motif_analyzer.py`: Position + sequence variation + cross-genome effects
  - `validate_consolidation.py`: Standalone validator (duplicated within consolidator CLI)
- `pipelines/`
  - `streme_pipeline.py`: Pipeline entry (prepare, consolidate, validate, analyze, full)
  - `Pipeline_promotor_discovery.yaml`: Genome-generic configuration reference
- `scripts/`
  - `remote_streme_all_lines.sh`: SLURM array job — one genome per task via `prepare`
  - `consolidate_streme_results.sh`: Example shell consolidation helper
- `docs/`: Additional documentation
  - `METHODS.md`: **Comprehensive Materials and Methods** for publication


## Inputs and expected files

**To start from genomes:** one or many `(genome, annotation)` pairs. Supply them on the command line (single genome) or in a tab-separated manifest:

```text
genome       fasta                          annotation                     expression
P_virgatus   /data/Penstemon_virgatus.fa    /data/Penstemon_virgatus.gff3  /data/virgatus_expr.tsv
P_barbatus   /data/Penstemon_barbatus.fa    /data/Penstemon_barbatus.gff3
```

Annotations may be GFF3 or GTF. The `expression` column is optional, as are `chromosomes` / `contig_pattern` columns for per-genome sequence filtering (see *Restricting to chromosomes* under Stage 0 below).

**To start from existing STREME results:** a directory containing one subfolder per genome, e.g.:

```text
streme_results/
├── streme_P_virgatus/sites.tsv
├── streme_P_barbatus/sites.tsv
└── streme_P_strictus/sites.tsv
```

Each `sites.tsv` should include (STREME defaults) columns similar to:

- `motif_ID`, `seq_ID`, `site_Start`, `site_End`, `site_Strand`, `site_Score`, `site_Sequence`

Genome names are derived from the `streme_<genome>` directory names; pass `--name-pattern` if your naming differs.

Expression files (tab-separated, two formats supported; the legacy `Line`/`line` headers are still accepted):

**Long format (Gene | Genome | Expression):**
```text
Gene       Genome      Expression
AT1G01010  P_barbatus  0.37
AT1G01010  P_strictus  -0.12
...
```

**Wide format (gene | LRTadd | <genome_1> | <genome_2> | ...):**
```text
gene       LRTadd  P_barbatus  P_strictus  P_virgatus
AT1G01010  2.45    0.37       -0.12        0.00
...
```

For relative analysis the reference genome serves as the baseline (expression = 0). It is auto-detected from an all-zero column, or set it explicitly with `--reference-genome`.


## Stage 0: Genome preparation (genomes → STREME results)

Turn genome assemblies + annotations into `streme_<genome>/` results. Multiple genomes run in parallel (`--jobs`); each genome's heavy steps use `--threads`.

```bash
# Many genomes from a manifest
python pipelines/streme_pipeline.py prepare \
  --manifest genomes.tsv --jobs 4 --threads 8 \
  --upstream 1000 --mask repeatmasker --species "Penstemon" \
  --output prepared/

# A single genome
python pipelines/streme_pipeline.py prepare \
  --genome P_virgatus --fasta P_virgatus.fa --annotation P_virgatus.gff3 \
  --upstream 1000 --threads 10 --output prepared/
```

Each genome flows through: **extract promoters** (1 kb upstream of TSS, strand-aware, pure Python) → **mask** (RepeatMasker or dust) → **background** (Markov model) → **STREME**. Stop early with `--mask none`, `--no-background`, or `--no-streme` (e.g. to only extract promoters, which needs no external tools). The individual steps are also available standalone via `cli_tools/genome_prep.py extract-promoters | mask | background | run-streme`.

The external tools are looked up on `PATH` by default. If they live at a module path (common on HPC), point at them explicitly with `--masker-path` (RepeatMasker/dust), `--fasta-get-markov-path`, and `--streme-path` — accepted on `prepare`, `full`, and the standalone `genome_prep.py` subcommands. A missing/unrunnable executable fails with a clear message.

**Custom RepeatMasker library (recommended for non-model organisms).** Dfam often ships only a tiny root partition and lacks plant clades; `--species "Penstemon"` will silently fall back to defaults (typically `homo sapiens`) and fail. Build a species-specific library once with RepeatModeler, then point `prepare` at it via `--mask-lib`:

```bash
# 1) One-time, per species (hours-to-days for a plant genome; runs BuildDatabase + RepeatModeler)
python cli_tools/genome_prep.py model-repeats penstemon_eatonii/PeChr.BYU.final.fa \
  --output-dir repeats/P_eatonii --name P_eatonii --threads 16
#   -> repeats/P_eatonii/P_eatonii-families.fa

# 2) Use it on every prepare run for that genome
streme-parser prepare --genome P_eatonii \
  --fasta penstemon_eatonii/PeChr.BYU.final.fa \
  --annotation penstemon_eatonii/penstemon_eatonii_final_nuclear_cp.gff3 \
  --upstream 2000 --contig-pattern '^PeChr' --threads 10 \
  --mask repeatmasker --mask-lib repeats/P_eatonii/P_eatonii-families.fa \
  --output prepared_penstemon_eatonii
```

For multi-genome runs, set the library per row in the manifest (`lib` column) instead — `lib` overrides `--mask-lib`, and `lib` takes precedence over `species` if both are set. `model-repeats` accepts `--repeatmodeler-path` / `--builddatabase-path` if those tools live outside `PATH`, and `--no-ltr-struct` to skip the (slow) LTR discovery stage.

**Restricting to chromosomes:** assemblies often mix assembled chromosomes with scaffolds, under varying names (`PeChr1…`, `Chr1`, `chr01`, scaffolds like `JBCEGF010000009.1`). Limit extraction to the sequences you want with either `--chromosomes` (an explicit allowlist — a comma list or a file with one name per line) or `--contig-pattern` (a regex); a feature is kept only if it passes both. Because naming differs per genome, these can also be set **per genome** in the manifest via `chromosomes` / `contig_pattern` columns, which override the run-wide flags for that row.

```bash
# Keep the 8 Penstemon eatonii chromosomes, drop the JBCEGF... scaffolds
python pipelines/streme_pipeline.py prepare \
  --genome P_eatonii --fasta PeChr.BYU.final.fa --annotation P_eatonii.gff3 \
  --contig-pattern '^PeChr' --output prepared/
# (equivalently: --chromosomes PeChr1,PeChr2,PeChr3,PeChr4,PeChr5,PeChr6,PeChr7,PeChr8)
```

Chromosome filtering happens at extraction time (when the annotation's sequence name is still known); if you start from pre-built STREME results it no longer applies, since those are keyed by gene only.

**FIMO scan of the STREME motifs.** STREME → FIMO is the canonical MEME-Suite chain: STREME discovers motifs *de novo*; FIMO then locates each motif's hits in your sequences at a controlled p- or q-value. Two ways to run it:

```bash
# 1) Inline — add FIMO to a prepare run (one fimo_<genome>/ per genome)
python pipelines/streme_pipeline.py prepare --manifest genomes.tsv --jobs 4 \
  --run-fimo --fimo-thresh 0.05 --fimo-qv-thresh --fimo-max-strand

# 2) Standalone — scan an existing prepared/ tree (useful for tuning thresholds
#    or re-scanning a different sequence set without rerunning STREME)
python pipelines/streme_pipeline.py scan prepared/ --jobs 4 \
  --thresh 0.05 --qv-thresh --max-strand
```

Defaults follow MEME-Suite guidance: for each genome `scan` picks `streme_<genome>/streme.txt` as the motif file, `<genome>_prep/<genome>_promoters.masked` as the sequence file, and `<genome>_prep/background.txt` as `--bgfile` (falling back to the motifs' embedded background if no background.txt). Override any of these with `--sequence` / `--bgfile`. Promoters are already ≤1 kb (matching FIMO's recommended scan length), so the default p-value 0.0001 yields a manageable false-positive rate; for safer reporting use `--qv-thresh` with `--thresh 0.05` (q-value ≤ 0.05).

The result is `prepared/streme_<genome>/` directories, ready for consolidation below.

## Quickstart (most direct path)

1. Consolidate all genomes’ STREME sites into a single table:

```bash
python cli_tools/streme_sites_consolidator.py consolidate prepared/ \
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

**Option A: Absolute analysis (per-genome analysis)**

Analyzes each genome independently. Expression data can be in long format (`Gene | Genome | Expression`) or wide format (one column per genome).

```bash
python cli_tools/motif_expression_analyzer.py \
  outputs/consolidated_streme_sites.tsv \
  path/to/expression.tsv \
  --output outputs/absolute_analysis_results
```

**Option B: Relative analysis (comparing to a reference-genome baseline)**

Generates features comparing each genome to the reference genome. Requires wide format expression data.

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
bin/streme-parser prepare --manifest genomes.tsv --jobs 4 --threads 8 --output prepared/
bin/streme-parser consolidate prepared/ --output outputs/
bin/streme-parser validate outputs/consolidated_streme_sites.tsv
bin/streme-parser analyze outputs/consolidated_streme_sites.tsv expression.tsv --type relative --output results/
bin/streme-parser analyze outputs/consolidated_streme_sites.tsv expression.tsv --type relative --reference-genome P_virgatus --output results/

# Full pipeline (consolidate → validate → analyze)
bin/streme-parser full prepared/ expression.tsv --output outputs/ --analysis-type relative
bin/streme-parser full prepared/ expression.tsv --output outputs/ --analysis-type relative --reference-genome P_virgatus

# Full pipeline starting from genomes (runs prepare first into the given dir)
bin/streme-parser full prepared/ expression.tsv --manifest genomes.tsv --jobs 4 --threads 8 --analysis-type relative
```

**Pipeline commands available:**
- `prepare`: Genomes → promoters → mask → background → STREME (parallel across genomes; optionally also FIMO via `--run-fimo`)
- `scan`: FIMO-scan a `prepared/` tree — locate STREME-motif hits per genome at a controlled p/q-value
- `consolidate`: Consolidate STREME motifs across genomes
- `validate`: Validate motif consolidation quality
- `analyze`: Run motif-expression analysis (absolute or relative)
- `full`: Complete pipeline (optionally prepare first, then consolidate → validate → analyze)

**Expression analysis types:**
- `--type absolute`: Per-genome analysis (each genome analyzed independently)
- `--type relative`: Comparative analysis (each genome vs reference baseline) - **recommended**
- `--reference-genome`: Set reference genome for relative analysis (`--reference-line` still accepted)

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
  `genome`, `gene_id`, `start_pos`, `end_pos`, `strand`, `score`, `sequence`,
  `relative_position`, `total_motifs_in_gene`, `relative_position_fraction`, `cluster_size`, `merged_count`, `length`
- Summary file reports counts per genome and top motifs

From features (`*_features.csv`):

- Row = one gene×genome combination; columns include per-motif presence, counts, position stats, and sequence variation stats
- Companion `*_feature_descriptions.txt` and `*_summary.txt` describe feature definitions and matrix stats

From motif→expression analysis:

- Cross-validated model performance across linear/regularized/tree models
- Feature importance (per-model and aggregate heatmap)
- Predictions and residuals per gene×genome (`predictions.tsv`)

From comprehensive analysis:

- Layer-wise R² contributions (presence, position, variation, cross-genome)
- Combined model performance and top features


## HPC workflow (SLURM example)

`scripts/remote_streme_all_lines.sh` is a SLURM array job that prepares one genome per task (extract → mask → background → STREME) from a manifest. Set the array size to the number of genome rows and submit:

```bash
sbatch --array=1-<N_genomes> scripts/remote_streme_all_lines.sh
```

After all tasks complete, consolidate the `streme_<genome>/` directories:

```bash
python cli_tools/streme_sites_consolidator.py consolidate prepared/ \
  --output outputs/consolidated_streme_sites
```

## Methodology and Statistical Approach

This pipeline implements two complementary approaches for motif-expression analysis:

### **Absolute Analysis**
- Models expression as a function of motif features within each genetic line
- Uses traditional cis-regulatory theory with additive motif effects
- Suitable for understanding line-specific regulatory programs

### **Relative Analysis** (Recommended)
- Models expression differences as a function of regulatory differences from a reference line
- Controls for genetic background and trans-acting factors
- Directly addresses the biological question: "How do regulatory changes explain expression changes?"

### **Statistical Framework**
- **Multiple algorithms**: Linear, Ridge, Lasso, Random Forest, Gradient Boosting
- **Cross-validation**: GroupKFold to prevent data leakage (genes never split between train/test)
- **Feature engineering**: Presence/absence, counts, positions, sequence composition
- **Performance metrics**: Cross-validated R², RMSE, feature importance

📖 **For complete Materials and Methods** (suitable for publication): See [`docs/METHODS.md`](docs/METHODS.md)

The methods document includes:
- Detailed algorithmic descriptions and parameter justifications
- Statistical model selection rationale with literature support
- Quality control procedures and validation approaches
- Assumptions, limitations, and recommended extensions
- Complete reference list for methodological approaches


## Tips and troubleshooting

- Use Python 3.8+ and install dependencies from `requirements.txt`.
- The consolidator tolerates missing/NaN fields in STREME `sites.tsv` rows and skips incomplete lines.
- If `pandas` isn’t available, output writing falls back to a manual TSV writer.
- Overlap merging is enabled by default in the consolidator to reduce sliding-window artifacts.
- For the regression analyzer, ensure your feature file delimiter matches its expectation (TSV). If you used the feature generator’s CSV output, convert to TSV (see Quickstart).
- Expression file (long format) must have headers: `Gene`, `Genome`, `Expression` (legacy `Line` still accepted).


## Citation and license

If you use this in a publication, please cite STREME (MEME Suite) for motif discovery. Project license and citation information for this toolkit are TBD.


## Acknowledgements

Built around MEME Suite STREME outputs. Includes utility scripts for an internal regulatory analysis workflow.
