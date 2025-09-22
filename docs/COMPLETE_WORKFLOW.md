# Complete Motif-Expression Analysis Workflow

This guide shows you exactly how to use your STREME parser to analyze motif effects on gene expression.

## 📊 The Analysis Question

**"Which regulatory motifs in gene promoters explain expression differences across genetic lines?"**

## 🔄 Complete Workflow

### Step 1: Consolidate STREME Results
```bash
# Consolidate motifs from all your STREME runs
./bin/streme-parser consolidate /path/to/streme/results/ --output analysis/
```
**Output:** `analysis/consolidated_streme_sites.tsv`

### Step 2: Extract Binary Motif Features
```bash
# Create binary presence/absence features for each gene-line combination
./bin/streme-parser extract-features analysis/consolidated_streme_sites.tsv \
    --simple \
    --output-prefix analysis/motif_features
```
**Output:** `analysis/motif_features_matrix.tsv`

### Step 3: Prepare Your Expression Data

Create a file called `expression_data.tsv` with this format:
```
Gene	Line	Expression
AT1G01010	IM502	2.3
AT1G01010	IM664	-1.2
AT1G01010	IM1034	0.8
AT1G01020	IM502	-0.5
...
```

**Key points:**
- Expression = log2(fold change) relative to IM767
- IM767 entries are not needed (they're always 0)
- Gene names must match those in your STREME analysis

### Step 4: Run Motif-Expression Analysis
```bash
# Analyze which motifs predict expression differences
./bin/streme-parser analyze-expression \
    analysis/motif_features_matrix.tsv \
    expression_data.tsv \
    --output analysis/regression_results/
```

### Step 5: Interpret Results

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
./bin/streme-parser extract-features analysis/consolidated_streme_sites.tsv \
    --output-prefix analysis/detailed_features

./bin/streme-parser analyze-expression \
    analysis/detailed_features_matrix.tsv \
    expression_data.tsv \
    --detailed \
    --output analysis/detailed_results/
```

### Filter to Most Variable Motifs
```bash
# Focus on top 100 most variable motifs
./bin/streme-parser analyze-expression \
    analysis/motif_features_matrix.tsv \
    expression_data.tsv \
    --top-motifs 100 \
    --output analysis/top100_results/
```

## 📝 What to Report to Your Advisor

1. **"We identified X motifs that explain Y% of expression variance across lines"**
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
- Ensure line names are consistent (IM502 vs IM_502)

### "Low R² values"
- Try `--detailed` mode for richer features
- Check if you have enough genes with variable expression
- Consider using `--top-motifs` to focus on most relevant motifs

### "Models failing"
- Check for missing values in expression data
- Ensure sufficient sample size (>100 gene-line combinations recommended)

**This workflow transforms your STREME results into actionable regulatory insights! 🧬→📊**