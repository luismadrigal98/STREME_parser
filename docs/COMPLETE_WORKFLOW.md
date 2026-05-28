# Complete Motif-Expression Analysis Workflow

This guide shows you exactly how to use your STREME parser to analyze motif effects on gene expression. It is genome-generic: the same workflow applies to several genomes of a genus (e.g. *Penstemon virgatus*, *P. barbatus*) or to multiple ecotypes/inbred lines of one species (e.g. *Mimulus guttatus* IM lines).

> **Terminology:** the canonical term is **genome**. Earlier versions used "line"; the legacy `line`/`Line` column names and the `--reference-line` flag are still accepted, so existing Mimulus data needs no migration.

## 📊 The Analysis Question

**"Which regulatory motifs in gene promoters explain expression differences across genomes?"**

## 🔄 Complete Workflow

### Step 1: Prepare Genomes (Promoters → STREME)
```bash
# From genome assemblies + annotations to STREME results, one dir per genome.
# Per genome: extract promoters (1 kb upstream of TSS) → mask → background → STREME.
python3 pipelines/streme_pipeline.py prepare \
    --manifest genomes.tsv --jobs 4 --threads 8 --upstream 1000 --output prepared/
```
**Output:** one `prepared/streme_<genome>/` directory per genome (e.g. `prepared/streme_P_virgatus/`).

### Step 2: Consolidate STREME Results
```bash
# Consolidate motifs from all your STREME runs
./bin/streme-parser consolidate prepared/ --output analysis/
```
**Output:** `analysis/consolidated_streme_sites.tsv`

### Step 3: Extract Binary Motif Features
```bash
# Create binary presence/absence features for each gene-genome combination
python3 cli_tools/streme_sites_consolidator.py features \
    analysis/consolidated_streme_sites.tsv \
    --simple \
    --output-prefix analysis/motif_features
```
**Output:** `analysis/motif_features_matrix.tsv`

### Step 4: Prepare Your Expression Data

Create a file called `expression_data.tsv` with this format (the `Genome` column accepts the legacy `Line` name as an alias):
```
Gene	Genome	Expression
AT1G01010	P_virgatus	2.3
AT1G01010	P_barbatus	-1.2
AT1G01010	IM767	0.8
AT1G01020	P_virgatus	-0.5
...
```

**Key points:**
- Expression = log2(fold change) relative to the reference genome (default: IM767)
- Reference-genome entries are not needed (they're always 0)
- Gene names must match those in your STREME analysis

### Step 5: Run Motif-Expression Analysis
```bash
# Analyze which motifs predict expression differences.
# Use --type absolute|relative and --reference-genome (alias --reference-line).
./bin/streme-parser analyze \
    analysis/consolidated_streme_sites.tsv \
    expression_data.tsv \
    --type relative --reference-genome IM767 \
    --output analysis/regression_results/
```

### Step 6: Interpret Results

The analysis creates:
- **`analysis_report.md`** - Main findings and interpretation
- **`motif_importance.tsv`** - Which motifs matter most
- **`predictions.tsv`** - Actual vs predicted expression for each gene
- **Visualizations** - Model performance and motif importance plots

## 📈 What the Analysis Tells You

### Example Results Interpretation:

**If R² = 0.65:**
🎉 "Motif patterns explain 65% of expression variance! Strong regulatory signal."

**If R² = 0.30:**
📈 "Motif patterns explain 30% of expression variance. Moderate regulatory effects."

**If R² = 0.05:**
⚠️ "Weak predictive power. Consider other factors (enhancers, chromatin, etc.)."

### Important Motifs Analysis:
The tool will show you which specific motifs are most predictive:
```
Top Motifs:
- motif_123_GATAAG: Importance 0.156 (strong predictor)
- motif_087_TTTCCC: Importance 0.142 (strong predictor)
- motif_234_AAATTT: Importance 0.098 (moderate predictor)
```

## 🧬 Biological Interpretation

### High Importance Motifs → Candidate Regulatory Elements
1. **Look up motif sequences** in transcription factor databases
2. **Check if known TF binding sites** match your important motifs
3. **Correlate with known biology** - do these TFs make sense for your traits?

### Gene-Level Predictions
Use `predictions.tsv` to find:
- **Well-predicted genes**: Strong motif-expression relationships
- **Poorly predicted genes**: May have other regulatory mechanisms
- **Outliers**: Interesting cases for follow-up

## 🔍 Advanced Analysis Options

### Focus on Specific Gene Sets
```bash
# Analyze only stress-response genes
grep "stress\|heat\|cold" gene_list.txt > stress_genes.txt
# Filter your expression data to these genes first
```

### Use Detailed Motif Features
```bash
# Include position, count, and sequence features (not just presence/absence)
python3 cli_tools/streme_sites_consolidator.py features \
    analysis/consolidated_streme_sites.tsv \
    --output-prefix analysis/detailed_features

./bin/streme-parser analyze \
    analysis/consolidated_streme_sites.tsv \
    expression_data.tsv \
    --type relative \
    --output analysis/detailed_results/
```

### Filter to Most Variable Motifs
```bash
# Focus on top 100 most variable motifs (feature extraction step)
python3 cli_tools/streme_sites_consolidator.py features \
    analysis/consolidated_streme_sites.tsv \
    --top-motifs 100 \
    --output-prefix analysis/top100_features
```

## 📝 What to Report to Your Advisor

1. **"We identified X motifs that explain Y% of expression variance across genomes"**
2. **"The top 10 predictive motifs are: [list with sequences]"**
3. **"These motifs match known binding sites for: [TF families]"**
4. **"Genes with strong motif-expression relationships include: [examples]"**
5. **"This approach successfully links promoter sequence to expression phenotype"**

## 🚀 Next Steps for Publication

1. **Validate predictions** - pick 5-10 highly predicted genes and test experimentally
2. **Gene ontology analysis** - are well-predicted genes in specific pathways?
3. **Motif refinement** - can you improve motifs using expression as feedback?
4. **Cross-validation** - test predictions on independent expression data
5. **Integration** - combine with other 'omics data (chromatin, metabolomics)

## 🛠 Troubleshooting

### "No overlapping data found"
- Check gene name formats match between motif and expression files
- Ensure genome names are consistent (P_virgatus vs P-virgatus)

### "Low R² values"
- Try detailed features (omit `--simple`) for richer inputs
- Check if you have enough genes with variable expression
- Consider using `--top-motifs` to focus on most relevant motifs

### "Models failing"
- Check for missing values in expression data
- Ensure sufficient sample size (>100 gene-genome combinations recommended)

**This workflow transforms your STREME results into actionable regulatory insights! 🧬→📊**