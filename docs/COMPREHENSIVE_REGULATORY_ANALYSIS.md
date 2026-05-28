# Multi-Layer Regulatory Analysis: Position, Variation & Expression

This guide shows you how to analyze **three layers of regulatory effects** simultaneously. It is genome-generic: the same analysis applies to several genomes of a genus (e.g. *Penstemon virgatus*, *P. barbatus*) or to multiple ecotypes/inbred lines of one species (e.g. *Mimulus guttatus* IM lines).

> **Terminology:** the canonical term is **genome**. Earlier versions used "line"; the legacy `line`/`Line` column names and the `--reference-line` flag are still accepted. The comprehensive analyzer's `cross_line` layer name is kept as a literal flag value for backward compatibility, but it measures cross-genome variation.

## 🎯 The Three Key Questions

1. **Motif Presence Effects**: Does having this motif affect expression?
2. **Position Bias Effects**: Does WHERE the motif is located matter for expression?
3. **Sequence Variation Effects**: Does HOW the motif varies across genomes affect expression?

## 🧬 Biological Hypothesis

**"Gene expression differences across genomes result from a combination of:**
- **Which motifs are present** (basic regulatory elements)
- **Where motifs are positioned** (distance from TSS matters)  
- **How motifs vary in sequence** (different variants have different strengths)"

## 📊 Complete Analysis Workflow

### Step 1: Prepare, Consolidate and Validate
```bash
# Prepare genomes (extract promoters → mask → background → STREME), one dir per genome
python3 pipelines/streme_pipeline.py prepare --manifest genomes.tsv --output prepared/

# Get your motif sites with exact positions and sequences, then validate
./bin/streme-parser consolidate prepared/ --output analysis/
./bin/streme-parser validate analysis/consolidated_streme_sites.tsv
```

### Step 2: Comprehensive Multi-Layer Analysis
```bash
# Analyze all regulatory layers together (run the analyzer directly)
python3 cli_tools/comprehensive_motif_analyzer.py \
    analysis/consolidated_streme_sites.tsv \
    expression_data.tsv \
    --output comprehensive_results/ \
    --layers presence position variation cross_line \
    --interactions \
    --top-motifs 100
```

### Step 3: Layer-by-Layer Comparison
```bash
# Test each layer individually to see their relative contributions
python3 cli_tools/comprehensive_motif_analyzer.py \
    analysis/consolidated_streme_sites.tsv \
    expression_data.tsv \
    --output layer_comparison/ \
    --layers presence \
    --top-motifs 100

python3 cli_tools/comprehensive_motif_analyzer.py \
    analysis/consolidated_streme_sites.tsv \
    expression_data.tsv \
    --output layer_comparison/ \
    --layers position \
    --top-motifs 100

python3 cli_tools/comprehensive_motif_analyzer.py \
    analysis/consolidated_streme_sites.tsv \
    expression_data.tsv \
    --output layer_comparison/ \
    --layers variation \
    --top-motifs 100
```

## 🔬 What Each Layer Captures

### Layer 1: Motif Presence
**Features Generated:**
- `motif_123_present`: Binary (0/1) - is motif present?
- `motif_123_count`: Integer - how many copies?

**Question:** "Do genes with motif X have different expression than genes without it?"

### Layer 2: Position Bias  
**Features Generated:**
- `motif_123_proximal_density`: Hits per 100bp in 0-500bp from TSS
- `motif_123_core_density`: Hits per 100bp in 500-1000bp from TSS  
- `motif_123_distal_density`: Hits per 100bp in 1000-2000bp from TSS
- `motif_123_position_weighted_score`: Exponentially weighted by distance from TSS
- `motif_123_spacing_regularity`: How evenly spaced are multiple hits
- `motif_123_closest_hit_distance`: Distance of nearest hit to TSS

**Question:** "Does motif X have stronger effects when closer to the TSS?"

### Layer 3: Sequence Variation
**Features Generated:**
- `motif_123_consensus_similarity_mean`: How similar are hits to consensus
- `motif_123_sequence_diversity`: How many different sequences per gene
- `motif_123_gc_content_mean`: Average GC content of motif hits
- `motif_123_unique_sequences`: Number of distinct sequences

**Question:** "Do different sequence variants of motif X have different regulatory strengths?"

### Layer 4: Cross-Genome Variation
**Features Generated** (literal feature names retain the legacy `line` wording for backward compatibility, but measure cross-genome variation):
- `motif_123_cross_line_diversity`: How variable is motif across all genomes
- `motif_123_line_has_unique_variant`: Does this genome have a unique sequence variant
- `motif_123_dominant_sequence_diversity`: How different are the most common variants

**Question:** "Do genome-specific sequence variants explain genome-specific expression differences?"

## 📈 Interpreting Results

### Example Output Interpretation:

```
## Regulatory Layer Contributions
- **Presence**: R² = 0.15 ± 0.03 (200 features)
- **Position**: R² = 0.35 ± 0.05 (600 features)  
- **Variation**: R² = 0.28 ± 0.04 (800 features)
- **Cross-genome**: R² = 0.22 ± 0.03 (300 features)

## Comprehensive Model Performance
**Combined Model**: R² = 0.67 ± 0.04

## Top 10 Most Important Features
- **motif_087_position_weighted_score** (Position): 0.0234
- **motif_123_consensus_similarity_mean** (Variation): 0.0198  
- **motif_045_proximal_density** (Position): 0.0187
- **motif_234_line_has_unique_variant** (Cross-genome): 0.0156
```

### Key Insights:
1. **Position matters most** (R² = 0.35) - WHERE motifs are located is crucial
2. **Sequence variation is important** (R² = 0.28) - HOW motifs vary affects function  
3. **Combined model is powerful** (R² = 0.67) - all layers work together
4. **Specific findings**: Motif 087 position near TSS is most predictive

## 🧬 Biological Implications

### High Position Importance → "Core Promoter Architecture Matters"
- Genes are sensitive to exact motif positioning
- TSS-proximal motifs have stronger effects
- Regulatory grammar depends on spatial organization

### High Variation Importance → "Sequence Variants Have Different Activities"  
- Not all motif hits are equal
- Sequence changes within motifs alter binding affinity
- Genome-specific variants may explain expression differences

### High Cross-Genome Importance → "Genetic Background Effects"
- Same motif behaves differently in different genomes
- Epistatic interactions between motifs and genetic background
- Genome-specific regulatory evolution

## 🔍 Advanced Analysis Options

### Focus on Specific Regulatory Questions:

```bash
# Question: "Is position bias the main driver?"
python3 cli_tools/comprehensive_motif_analyzer.py data.tsv expr.tsv --layers position

# Question: "Do sequence variants explain genome differences?"  
python3 cli_tools/comprehensive_motif_analyzer.py data.tsv expr.tsv --layers variation cross_line

# Question: "Are there motif-motif interactions?"
python3 cli_tools/comprehensive_motif_analyzer.py data.tsv expr.tsv --interactions
```

### Gene Set Specific Analysis:
```bash
# Focus on stress-response genes
grep "stress\|heat\|drought" gene_annotations.txt > stress_genes.txt
# Filter your expression data to these genes first
```

## 📊 What to Report

### For Your Advisor:
1. **"We developed a multi-layer regulatory model that explains X% of expression variance"**
2. **"Position effects account for Y% - proximity to TSS is crucial"**  
3. **"Sequence variation accounts for Z% - motif variants have different activities"**
4. **"The top predictive features are: [specific motifs and their properties]"**
5. **"This approach successfully decomposes regulatory effects into interpretable components"**

### For Publication:
1. **Model performance comparison** across layers
2. **Feature importance analysis** showing which motifs/positions/variants matter
3. **Specific examples** of genes well-predicted by the model
4. **Biological validation** of top predictive motifs
5. **Cross-validation** showing robustness across gene sets

## 🎯 Expected Outcomes

### Strong Model (R² > 0.5):
✅ **"Regulatory sequence successfully predicts expression"**
- Clear biological signal in promoter architecture
- Publishable regulatory model
- Candidates for experimental validation

### Moderate Model (R² = 0.2-0.5):
📈 **"Partial regulatory explanation - other factors involved"**
- Promoter contributes but other mechanisms important
- Consider chromatin, enhancers, post-transcriptional regulation
- Still valuable for understanding promoter contributions

### Weak Model (R² < 0.2):
🔍 **"Limited promoter effects - look elsewhere"**
- Expression differences mainly from other sources
- Consider longer sequences, chromatin state, trans effects
- May indicate your trait is not primarily promoter-driven

This comprehensive approach gives you the complete regulatory picture! 🧬📊