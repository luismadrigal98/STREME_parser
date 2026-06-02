# STREME Analysis Pipeline

A comprehensive toolkit for analyzing STREME motif discovery results across multiple genomes. It is genome-generic: use it for several genomes of a genus (e.g. *Penstemon virgatus*, *P. barbatus*, *P. strictus*) or for multiple ecotypes/inbred lines of one species (e.g. *Mimulus guttatus* IM lines).

> **Terminology:** the canonical term is **genome**. Earlier versions used "line"; the legacy `line`/`Line` column names and the `--reference-line` flag are still accepted everywhere, so existing Mimulus data needs no migration.

## Quick Start

The pipeline provides a unified CLI tool that can run different analysis subprograms:

```bash
# Prepare genome(s): extract promoters, mask, build background, run STREME
python3 pipelines/streme_pipeline.py prepare \
  --manifest genomes.tsv --jobs 4 --threads 8 --upstream 1000 --output prepared/

# Consolidate motifs across genomes
./bin/streme-parser consolidate prepared/ --output my_analysis

# Extract features for machine learning (consolidator CLI)
python3 cli_tools/streme_sites_consolidator.py features my_analysis/consolidated_streme_sites.tsv --simple

# Or run any pipeline step using Python directly
python3 pipelines/streme_pipeline.py consolidate prepared/ --output my_analysis
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

### `prepare` - Genome Preparation (Promoters → STREME)

Goes from genome assemblies plus annotations to STREME results. Per genome the chain is: extract promoters (N bp upstream of the TSS, default 1000, strand-aware, pure Python) → mask (RepeatMasker or dust) → background (Markov model) → STREME, producing a `prepared/streme_<genome>/` directory for each genome.

```bash
# Many genomes from a manifest, run in parallel (--jobs genomes at a time)
python3 pipelines/streme_pipeline.py prepare \
  --manifest genomes.tsv --jobs 4 --threads 8 --upstream 1000 \
  --mask repeatmasker --species "Penstemon" --output prepared/

# A single genome (parallelise its heavy steps with --threads)
python3 pipelines/streme_pipeline.py prepare \
  --genome P_virgatus --fasta P_virgatus.fa --annotation P_virgatus.gff3 \
  --upstream 1000 --threads 10 --output prepared/
```

The `--manifest` is a tab-separated file with a header row and the columns `genome`, `fasta`, `annotation`, and optional `expression`, `chromosomes`, and `contig_pattern`:

```
genome      fasta                 annotation             expression               contig_pattern
P_virgatus  /data/virgatus.fa     /data/virgatus.gff3    /data/virgatus_expr.tsv  ^Chr
P_eatonii   /data/eatonii.fa      /data/eatonii.gff3     /data/eatonii_expr.tsv   ^PeChr
P_barbatus  /data/barbatus.fa     /data/barbatus.gff3
```

**Restricting to chromosomes:** assemblies often mix chromosomes with scaffolds under different names (`PeChr1…`, `Chr1`, `JBCEGF010000009.1`). Limit extraction with `--chromosomes` (a comma list or a file, one name per line) and/or `--contig-pattern` (a regex); a feature is kept only if it passes both. Because naming varies per genome, set these per genome via the manifest `chromosomes` / `contig_pattern` columns, which override the run-wide flags for that row. Filtering applies at extraction time (it relies on the annotation's sequence name).

Standalone steps are also available via `cli_tools/genome_prep.py`:

```bash
# Extract promoters only (pure Python, no external tools needed)
python3 cli_tools/genome_prep.py extract-promoters genome.fa annot.gff3 \
  -o promoters.fa --upstream 1000 --feature-type gene --avoid-overlap \
  --contig-pattern '^PeChr'   # keep chromosomes, drop scaffolds

# Plus: mask, background, run-streme subcommands
```

The `mask`, `background`, `run-streme`, and `run-fimo` steps auto-detect the external tool they need on `PATH` and error clearly if it is missing. If a tool lives at a module path instead (common on HPC), point at it with `--masker-path` (RepeatMasker/dust), `--fasta-get-markov-path`, `--streme-path`, or `--fimo-path` — available on `prepare`, `full`, `scan`, and the standalone subcommands.

**Custom RepeatMasker library (non-model organisms).** Dfam's bundled root partition often lacks plant clades, so `--species` will fail. Build a species-specific library once with RepeatModeler, then point `prepare` at it via `--mask-lib`:

```bash
# 1) One-time per species (hours-to-days)
python3 cli_tools/genome_prep.py model-repeats penstemon_eatonii/PeChr.BYU.final.fa \
  --output-dir repeats/P_eatonii --name P_eatonii --threads 16
#   -> repeats/P_eatonii/P_eatonii-families.fa

# 2) Use on every prepare run
streme-parser prepare --genome P_eatonii --fasta ... --annotation ... \
  --mask repeatmasker --mask-lib repeats/P_eatonii/P_eatonii-families.fa \
  --threads 10 --output prepared/
```

For multi-genome runs set the library per row via the manifest `lib` column. `lib` overrides `--mask-lib`, and `lib` takes precedence over `species` when both are present for a genome.

**Tandem repeats / SSRs.** RepeatMasker libraries cover TEs but not simple-sequence repeats; the leftover `(AT)n`, `(GAA)n`, etc. tracts in plant promoters routinely dominate STREME output as spurious AT-rich motifs. Chain Tandem Repeats Finder (`bioconda::trf`) after the primary masker with `--extra-mask trf` (or `dust`, or `both` to chain dust then TRF). The pipeline writes `<genome>_promoters.masked` → `.masked.dust` → `.masked.trf`, and STREME reads the final file. Use `--trf-path` if `trf` isn't on PATH.

### `scan` - FIMO motif scan (STREME → FIMO)

STREME discovers motifs *de novo*; FIMO then locates each motif's hits in your sequences at a controlled p- or q-value. Run inline as part of `prepare`, or as a standalone step against an existing `prepared/` tree:

```bash
# Inline: add --run-fimo to a prepare run
python3 pipelines/streme_pipeline.py prepare --manifest genomes.tsv --jobs 4 \
  --run-fimo --fimo-thresh 0.05 --fimo-qv-thresh --fimo-max-strand

# Standalone: scan a previously prepared tree (handy for tuning thresholds)
python3 pipelines/streme_pipeline.py scan prepared/ --jobs 4 \
  --thresh 0.05 --qv-thresh --max-strand
```

For each `streme_<genome>/` the scan auto-picks `streme.txt` as the motif file, `<genome>_prep/<genome>_promoters.masked` as the sequence file, and `<genome>_prep/background.txt` as `--bgfile` (falling back to the motifs' embedded background). Override any with `--sequence`/`--bgfile`. Output goes to `fimo_<genome>/`.

### `annotate` - TOMTOM motif vs reference TF database

Matches each genome's STREME motifs against a reference TF PWM database so the de novo motifs get putative TF assignments. Inline via `prepare --run-tomtom --tomtom-db <db>` or standalone:

```bash
python3 pipelines/streme_pipeline.py annotate prepared/ \
  --target-db ~/dbs/Ath_TF_binding_motifs.meme --jobs 4 --thresh 0.1
```

Output: `tomtom_<genome>/` next to each `streme_<genome>/`. Default `--thresh 0.1` is a q-value cutoff; switch to E-value with `--evalue`. Common downloadable databases (MEME format): **PlantTFDB 5.0** (Arabidopsis), **JASPAR plants**, **CIS-BP**. Use `--tomtom-path` if `tomtom` isn't on PATH.

### `network` - exploratory analysis without expression data

When you don't have expression data — or want to look at the regulatory landscape independently of it — `network` runs two complementary analyses on a consolidated TSV:

- **Motif co-occurrence network**: motif pairs scored by Jaccard / lift / Fisher one-sided p (BH-corrected); pairs passing thresholds become edges; connected components → motif modules (combinatorial regulation hypotheses). Exported as **GraphML** for Cytoscape / Gephi.
- **Gene clustering by motif profile**: hierarchical Jaccard / average-linkage clustering of genes; per-cluster motif **fingerprints** (Fisher enrichment vs background) summarise what defines each cluster.
- **Positional summary**: per-motif distance-from-promoter-start and relative-position distributions.

```bash
streme-parser consolidate prepared_penstemon_eatonii/ --output outputs/consolidated_P_eatonii
streme-parser validate  outputs/consolidated_P_eatonii.tsv
streme-parser network   outputs/consolidated_P_eatonii.tsv --output network_results/ --n-clusters 30
```

Outputs: `motif_cooccurrence_edges.tsv`, `motif_cooccurrence_network.graphml`, `motif_modules.tsv`, `gene_clusters.tsv`, `gene_cluster_fingerprints.tsv`, `motif_positional_summary.tsv`, `network_summary.txt`. After TOMTOM annotation, the motif IDs in modules/fingerprints can be joined to TF names from `tomtom.tsv` for biological interpretation.

### `consolidate` - Consolidate STREME Sites

Processes STREME `sites.tsv` files from multiple genomes, consolidates similar motifs, and creates comprehensive regulatory maps.

#### Key Features:
- ✅ **Motif Consolidation**: Groups similar motifs using IUPAC-aware sequence comparison
- ✅ **Overlap Merging**: Removes redundant overlapping hits (solves poly-A duplication issues)
- ✅ **Repetitive Sequence Handling**: Proper handling of poly-A, poly-T tracts
- ✅ **Exact Coordinates**: Preserves genomic positions and sequences
- ✅ **Comprehensive Output**: 17 columns with detailed motif information

#### Usage Examples:

```bash
# Basic consolidation with overlap merging (recommended)
./bin/streme-parser consolidate prepared/ --output comprehensive_analysis

# Custom similarity threshold for stricter motif grouping
./bin/streme-parser consolidate prepared/ \
  --threshold 0.8 \
  --output strict_analysis
```

The consolidator CLI (`cli_tools/streme_sites_consolidator.py`) exposes additional filtering and merging options:

```bash
# Process specific genomes only ("--lines" is an accepted legacy alias)
python3 cli_tools/streme_sites_consolidator.py consolidate prepared/ \
  --genomes P_virgatus,P_barbatus \
  --output subset_analysis

# Restrict input directories by name pattern
python3 cli_tools/streme_sites_consolidator.py consolidate prepared/ \
  --name-pattern 'streme_*' \
  --output pattern_analysis

# Custom overlap threshold for merging
python3 cli_tools/streme_sites_consolidator.py consolidate prepared/ \
  --overlap-threshold 0.7 \
  --output custom_overlap

# Disable overlap merging (for debugging)
python3 cli_tools/streme_sites_consolidator.py consolidate prepared/ \
  --no-merge-overlaps \
  --output debug_analysis

# Verbose output for troubleshooting
python3 cli_tools/streme_sites_consolidator.py consolidate prepared/ \
  --verbose \
  --output verbose_analysis
```

#### Input Directory Structure:
```
prepared/
├── streme_P_virgatus/sites.tsv
├── streme_P_barbatus/sites.tsv
├── streme_P_strictus/sites.tsv
└── streme_IM767/sites.tsv
```

#### Output Files:
- `{output_prefix}.tsv` - Consolidated motif sites with 17 columns
- `{output_prefix}_summary.txt` - Statistical summary of results

### `features` - Extract Regression Features

Feature extraction is provided by the consolidator CLI (`cli_tools/streme_sites_consolidator.py features`). It converts consolidated motif data into machine-learning-ready features for gene expression prediction.

#### Key Features:
- ✅ **Binary Mode**: Simple presence/absence features for each motif (--simple flag)
- ✅ **Detailed Mode**: Comprehensive features including positional, sequence variation, and density metrics
- ✅ **Expression Integration**: Optional expression data integration
- ✅ **Motif Filtering**: Filter by frequency or select top N motifs
- ✅ **Ready for ML**: Output format compatible with scikit-learn and other ML libraries

#### Usage Examples:

```bash
# Simple binary features (presence/absence only)
python3 cli_tools/streme_sites_consolidator.py features consolidated_streme_sites.tsv --simple

# Detailed features with expression data
python3 cli_tools/streme_sites_consolidator.py features consolidated_streme_sites.tsv \
  --expression expression_data.tsv \
  --top-motifs 100

# Filter motifs by minimum occurrence
python3 cli_tools/streme_sites_consolidator.py features consolidated_streme_sites.tsv \
  --min-sites 50 \
  --output-prefix filtered_features
```

#### Feature Types:
- **Simple Mode (--simple)**: Only binary presence/absence features
- **Detailed Mode**: Positional features, sequence variation, GC content, density metrics, spacing regularity

### Expression Analysis

After consolidation, relate motifs to expression with the pipeline `analyze` command. It takes `--type absolute|relative` and `--reference-genome` (alias `--reference-line`):

```bash
# Relative analysis vs a reference genome (default reference: IM767)
./bin/streme-parser analyze comprehensive_analysis.tsv expression_data.tsv \
  --type relative --reference-genome P_strictus --output results/
```

## Output Columns

The consolidated TSV file contains 17 columns with comprehensive motif information:

| Column | Description |
|--------|-------------|
| `consolidated_motif_id` | Unified motif identifier (MOTIF_001, etc.) |
| `original_motif_id` | Original STREME identifier (1-AAAAA, etc.) |
| `motif_consensus` | True consensus computed from observed sequences |
| `original_streme_consensus` | Original STREME consensus pattern (prefix removed) |
| `genome` | Genome identifier (legacy alias: `line`) |
| `gene_id` | Gene identifier |
| `start_pos`, `end_pos` | Genomic coordinates |
| `strand` | DNA strand (+ or -) |
| `score` | STREME confidence score |
| `sequence` | Actual sequence found |
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
MOTIF_001  IM767  MiIM7v11000019m.g  786  800  AAAATAAAAAAAAAG
MOTIF_001  IM767  MiIM7v11000019m.g  787  801  AAAAATAAAAAAAAA
MOTIF_001  IM767  MiIM7v11000019m.g  788  802  AAAAAATAAAAAAAA
```

**After** (merged):
```
MOTIF_001  IM767  MiIM7v11000019m.g  786  805  AAAATAAAAAAAAAG  merged_count: 5
```

### Improved Similarity Algorithm
- Better handling of repetitive sequences (poly-A, poly-T)
- Reduced length penalties for biological repeats
- Maintains precision for complex motifs

## Documentation

- **[`OUTPUT_COLUMNS_GUIDE.md`](OUTPUT_COLUMNS_GUIDE.md)** - Detailed column descriptions
- **[`METHODS.md`](METHODS.md)** - Materials and Methods (publication-oriented methodology)
- **[`../README.md`](../README.md)** - Project overview and complete usage examples

## Typical Results

For 4 genomes with ~1000 genes each:
- **Input**: ~100,000 raw STREME motif sites
- **After consolidation**: ~180 unique motif patterns
- **After overlap merging**: ~80,000 non-redundant sites
- **Processing time**: 2-5 minutes

## Troubleshooting

### Common Issues:

1. **Missing pandas**: `pip install pandas`
2. **Permission denied**: `chmod +x bin/streme-parser`
3. **No sites.tsv files found**: Check directory structure
4. **Empty output**: Use `--verbose` flag for debugging

### Getting Help:

```bash
./bin/streme-parser --help                    # General help
./bin/streme-parser consolidate --help        # Consolidation-specific help
```

## Contributing

The pipeline is designed to be modular and extensible. New analysis tools can be added as subcommands in `pipelines/streme_pipeline.py`.

## Example Workflow

```bash
# 1. Prepare genomes and consolidate all STREME results
python3 pipelines/streme_pipeline.py prepare --manifest genomes.tsv --output prepared/
./bin/streme-parser consolidate prepared/ \
  --output comprehensive_analysis \
  --verbose

# 2. Examine results
head comprehensive_analysis.tsv
cat comprehensive_analysis_summary.txt

# 3. Find genes with high regulatory complexity
awk -F'\t' '$13 > 10' comprehensive_analysis.tsv | cut -f6 | sort | uniq

# 4. Analyze motif distribution patterns
cut -f1 comprehensive_analysis.tsv | sort | uniq -c | sort -nr | head -20
```

This provides a **clean, unified interface** for all your motif analysis needs!
