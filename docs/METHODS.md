# Materials and Methods

## Motif Discovery and Consolidation

### STREME Motif Discovery
Motif discovery was performed using STREME (STochastic Regular Expression Motif Elicitation) version 5.4.1 from the MEME Suite (Bailey et al., 2015; Bailey, 2021). STREME was run independently on promoter sequences for each genetic line using default parameters with the following settings: minimum motif width of 6 bp, maximum motif width of 15 bp, and a maximum of 20 motifs per analysis. The algorithm was configured to identify both palindromic and non-palindromic motifs using a first-order Markov model for background sequence generation.

### Motif Consolidation Across Genetic Lines
To address the challenge of comparing motifs discovered independently across multiple genetic lines, we implemented a comprehensive motif consolidation pipeline that accounts for the biological variability in motif instances while maintaining statistical rigor.

#### IUPAC-Aware Sequence Similarity Calculation
Motif similarity was calculated using IUPAC (International Union of Pure and Applied Chemistry) nucleotide ambiguity codes to properly handle degenerate positions. For each pair of motifs, we computed alignment scores using the following approach:

1. **Dynamic Programming Alignment**: Motifs were aligned using the Needleman-Wunsch algorithm modified for IUPAC codes, where matches between compatible nucleotides (e.g., 'A' with 'R' [A or G]) received partial scores based on overlap probability.

2. **Length Penalty Application**: To prevent spurious matches between motifs of substantially different lengths, we applied a length penalty factor:
   ```
   Length_Penalty = 1 - |len(motif1) - len(motif2)| / max(len(motif1), len(motif2))
   ```

3. **Similarity Score Calculation**: The final similarity score was computed as:
   ```
   Similarity = (Alignment_Score / max_possible_score) × Length_Penalty
   ```

#### Clustering and Consensus Generation
Motifs with similarity scores ≥ 0.75 (determined through empirical validation to balance specificity and sensitivity) were clustered using single-linkage clustering. For each cluster:

1. **True Consensus Calculation**: Rather than relying solely on STREME's reported consensus sequences, we computed true consensus sequences from the actual binding site instances, providing more accurate representations of the discovered motifs.

2. **Cluster Quality Assessment**: Each cluster was evaluated for coherence using within-cluster similarity metrics and manual inspection of representative sequences.

### Site Overlap Resolution and Merging
Within each gene promoter, overlapping motif sites were resolved using a priority-based merging system:

1. **Overlap Detection**: Sites were considered overlapping if they shared ≥50% of their sequence positions.

2. **Priority Assignment**: In cases of overlap, sites were prioritized based on: (1) STREME significance scores, (2) motif cluster size (larger clusters indicating more widespread occurrence), and (3) sequence length.

3. **Coordinate Adjustment**: Merged sites retained the boundaries of the highest-priority contributing site, with updated significance scores reflecting the combined evidence.

## Expression Data Processing and Analysis

### Data Preprocessing
Gene expression data were processed to support both long format (Gene | Line | Expression) and wide format (gene | LRTadd | IM62 | IM155 | ... | IM767) inputs. For relative analysis, expression values represent log2 fold-changes relative to the designated reference line (default: IM767), ensuring that the reference line has expression values of zero across all genes.

### Feature Engineering
We developed a comprehensive feature engineering pipeline that captures multiple dimensions of regulatory variation:

#### Presence/Absence Features
- **Binary presence**: Motif present (1) or absent (0) in gene promoter
- **Count features**: Number of motif instances per gene
- **Density features**: Motif instances per kb of promoter sequence

#### Positional Features
- **Position mean**: Average distance from transcription start site (TSS)
- **Position variance**: Variability in motif positioning within promoter
- **Proximal density**: Density of motifs within 500bp of TSS (Lenhard et al., 2012)
- **Relative position fraction**: Motif position as fraction of total promoter length

#### Sequence Composition Features
- **GC content**: Average GC content of motif instances
- **Sequence diversity**: Shannon entropy of nucleotide composition within motif instances
- **Degeneracy index**: Measure of IUPAC ambiguity within motif consensus

### Statistical Modeling Framework

#### Model Selection Rationale
We employed an ensemble approach using five complementary machine learning algorithms, each capturing different aspects of motif-expression relationships:

1. **Linear Regression**: Baseline model assuming additive motif effects (Montgomery et al., 2010)
2. **Ridge Regression**: L2-regularized linear model preventing overfitting with correlated features (Hoerl & Kennard, 1970)
3. **Lasso Regression**: L1-regularized model for automatic feature selection (Tibshirani, 1996)
4. **Random Forest**: Non-parametric ensemble method capturing non-linear interactions (Breiman, 2001)
5. **Gradient Boosting**: Sequential ensemble method for complex pattern recognition (Friedman, 2001)

#### Cross-Validation Strategy
To prevent data leakage and ensure robust performance estimates, we implemented GroupKFold cross-validation with genes as grouping units (Varoquaux et al., 2017). This approach ensures that:
- The same gene never appears in both training and testing sets
- Model performance reflects ability to generalize to unseen genes
- Estimates are not inflated by within-gene correlation structure

#### Performance Metrics
Model performance was evaluated using multiple complementary metrics:
- **Cross-validated R²**: Proportion of variance explained in held-out data
- **Root Mean Square Error (RMSE)**: Absolute prediction accuracy
- **Feature importance**: Variable importance scores across models

### Absolute vs. Relative Analysis Frameworks

#### Absolute Analysis
The absolute analysis framework models expression as a direct function of motif features within each genetic line:

```
Expression_ij = β₀ + Σₖ(βₖ × Motif_Feature_ijk) + εᵢⱼ
```

Where:
- i = gene index
- j = genetic line index  
- k = motif feature index
- ε = residual error term

This approach is based on classical cis-regulatory theory where transcription factor binding sites additively contribute to gene expression (Bintu et al., 2005; Segal et al., 2008).

#### Relative Analysis (Novel Approach)
The relative analysis framework directly models expression differences as a function of regulatory feature differences from a reference line:

```
ΔExpression_ij = β₀ + Σₖ(βₖ × ΔMotif_Feature_ijk) + εᵢⱼ
```

Where:
- ΔExpression_ij = Expression_ij - Expression_i,reference
- ΔMotif_Feature_ijk = Motif_Feature_ijk - Motif_Feature_i,reference,k

**Theoretical Justification**: This approach is grounded in evolutionary biology and comparative genomics, where regulatory evolution is understood through changes relative to ancestral states (Wittkopp & Kalay, 2012). By using an isogenic reference line, we control for:
- Genetic background effects
- Trans-acting factors common across lines
- Technical batch effects in expression measurement
- Baseline chromatin accessibility differences

**Statistical Advantages**:
1. **Reduced confounding**: Reference-based normalization removes systematic biases
2. **Direct biological relevance**: Models the actual biological question (what changes cause expression differences?)
3. **Improved statistical power**: Focuses on variation of interest rather than absolute levels

### Feature Importance and Interpretation

#### Multi-Model Consensus
Feature importance was calculated using algorithm-specific methods:
- **Linear models**: Absolute regression coefficients
- **Tree-based models**: Gini impurity reduction or permutation importance
- **Consensus ranking**: Features ranked by average importance across all models

#### Statistical Significance
Rather than relying on traditional p-values (which can be misleading in high-dimensional settings; Efron, 2010), we emphasize:
- **Cross-validated performance**: Robust measure of predictive value
- **Consistency across models**: Features important in multiple algorithms
- **Effect size magnitude**: Biological significance of observed associations

### Quality Control and Validation

#### Model Validation
1. **Overfitting Detection**: Large gaps between training and cross-validation performance
2. **Residual Analysis**: Systematic patterns in prediction residuals
3. **Feature Stability**: Consistency of important features across cross-validation folds

#### Biological Validation Recommendations
1. **Literature Concordance**: Comparison with known transcription factor binding sites
2. **Functional Enrichment**: Gene Ontology analysis of genes with high motif content
3. **Experimental Validation**: ChIP-seq or reporter assay confirmation of predicted regulatory relationships

## Computational Implementation

### Software and Dependencies
The analysis pipeline was implemented in Python 3.8+ using the following packages:
- **pandas** (1.3+): Data manipulation and analysis
- **numpy** (1.20+): Numerical computing
- **scikit-learn** (1.0+): Machine learning algorithms
- **matplotlib/seaborn**: Data visualization
- **MEME Suite** (5.4.1): Motif discovery (Bailey et al., 2015)

### Computational Resources and Scalability
The pipeline is designed for computational efficiency:
- **Memory usage**: Linear scaling with number of genes and motifs
- **Parallel processing**: Cross-validation and ensemble methods utilize multiple cores
- **Batch processing**: Support for large-scale genomic datasets

### Reproducibility
All analyses include:
- **Fixed random seeds**: Ensuring reproducible machine learning results
- **Version control**: Complete parameter and software version tracking  
- **Intermediate file preservation**: Full audit trail of analysis steps

## Statistical Assumptions and Limitations

### Assumptions
1. **Independence**: Gene expression values are independent after accounting for genetic line
2. **Linearity**: Motif effects are approximately additive (relaxed in tree-based models)
3. **Stationarity**: Regulatory relationships are consistent across the analyzed conditions
4. **Completeness**: Analyzed motifs capture major regulatory variation

### Limitations
1. **Correlation vs. Causation**: Associations do not imply direct regulatory causation
2. **Context Independence**: Does not model tissue- or condition-specific regulatory interactions
3. **Single-layer Analysis**: Does not account for hierarchical gene regulatory networks
4. **Motif Discovery Bias**: Limited to motifs detectable by STREME algorithm

### Recommended Extensions
1. **Experimental validation** through ChIP-seq or functional genomics approaches
2. **Network analysis** incorporating known transcription factor interactions
3. **Multi-condition analysis** examining regulatory plasticity across environments
4. **Integration with chromatin accessibility** data (ATAC-seq or DNase-seq)

## References

Bailey, T.L. (2021). STREME: accurate and versatile sequence motif discovery. *Bioinformatics*, 37(18), 2834-2840.

Bailey, T.L., et al. (2015). The MEME Suite. *Nucleic Acids Research*, 43(W1), W39-W49.

Bintu, L., et al. (2005). Transcriptional regulation by the numbers: models. *Current Opinion in Genetics & Development*, 15(2), 116-124.

Breiman, L. (2001). Random forests. *Machine Learning*, 45(1), 5-32.

Efron, B. (2010). Large-scale inference: empirical Bayes methods for estimation, testing, and prediction. Cambridge University Press.

Friedman, J.H. (2001). Greedy function approximation: a gradient boosting machine. *Annals of Statistics*, 1189-1232.

Hoerl, A.E., & Kennard, R.W. (1970). Ridge regression: Biased estimation for nonorthogonal problems. *Technometrics*, 12(1), 55-67.

Lenhard, B., et al. (2012). Metazoan promoters: emerging characteristics and insights into transcriptional regulation. *Nature Reviews Genetics*, 13(4), 233-245.

Montgomery, S.B., et al. (2010). Transcriptome genetics using second generation sequencing in a Caucasian population. *Nature*, 464(7289), 773-777.

Segal, E., et al. (2008). Predicting expression patterns from regulatory sequence in Drosophila segmentation. *Nature*, 451(7178), 535-540.

Tibshirani, R. (1996). Regression shrinkage and selection via the lasso. *Journal of the Royal Statistical Society*, 58(1), 267-288.

Varoquaux, G., et al. (2017). Cross-validation failure: Small sample sizes lead to large error bars. *NeuroImage*, 180, 68-77.

Wittkopp, P.J., & Kalay, G. (2012). Cis-regulatory elements: molecular mechanisms and evolutionary processes underlying divergence. *Nature Reviews Genetics*, 13(1), 59-69.