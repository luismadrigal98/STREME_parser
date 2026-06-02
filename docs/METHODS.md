# Materials and Methods

The pipeline is genome-generic and operates on any set of genome assemblies with accompanying annotations. It has been applied both to multiple genomes of a genus (e.g. *Penstemon virgatus*, *P. barbatus*, *P. strictus*) and to multiple ecotypes/inbred lines of a single species (e.g. *Mimulus guttatus* IM lines). Throughout, the canonical unit of comparison is referred to as a "genome"; earlier versions of the toolkit used the term "line", which is retained as an accepted alias in column names and command-line flags so that legacy datasets require no migration.

## Promoter Extraction and Genome Preparation

### Promoter Sequence Extraction
For each genome, candidate promoter regions were extracted directly from the assembly using the gene annotation. By default, the 1 kb region immediately upstream of each gene's transcription start site (TSS) was retrieved. Extraction is strand-aware: for genes on the plus strand the window spans the bases 5' of the annotated start coordinate, whereas for genes on the minus strand the window is taken downstream of the end coordinate and reverse-complemented, so that all extracted sequences are oriented 5'→3' relative to the gene. Windows were clipped at contig boundaries to avoid running past the ends of scaffolds. The upstream window length (and an optional downstream extension past the TSS) is configurable, as is the annotation feature type used to define gene starts (gene by default). Optionally, promoter windows can be truncated short of the nearest adjacent gene to avoid overlap with neighbouring loci. Extraction can also be restricted to a chosen subset of sequences — typically the assembled chromosomes, excluding unplaced scaffolds — by an explicit sequence allowlist or a regular expression matched against sequence names; because naming conventions differ between assemblies (e.g. `PeChr1`–`PeChr8` versus scaffold accessions), this filter can be specified per genome. Promoter extraction is implemented in pure Python and requires no external tools.

### Repeat Masking and Background Modelling
Extracted promoter sequences were optionally masked for repetitive and low-complexity content using either RepeatMasker (Smit, Hubley & Green, 1996–2015) or the dust algorithm, and a Markov background model (first order by default) was estimated from the masked sequences for use by STREME. These steps detect the required external tool at run time and report a clear error if it is unavailable.

For non-model organisms whose lineages are absent from the bundled Dfam partitions, a species-specific repeat library was constructed once using RepeatModeler (Flynn et al., 2020) — wrapping the standard `BuildDatabase` + `RepeatModeler -LTRStruct` workflow — and supplied to RepeatMasker via its `-lib` option for all subsequent preparation runs of that genome. In multi-genome runs, each genome's library is provided independently through the per-genome `lib` column of the input manifest; a custom library takes precedence over any species-based Dfam lookup when both are configured.

### Genome-Parallel Preparation
The promoter extraction, masking, background-modelling, and STREME steps are orchestrated per genome by a single preparation command. Multiple genomes can be supplied through a manifest and processed in parallel, while the computationally intensive per-genome steps (masking and motif discovery) can themselves be multithreaded. The result is one STREME output directory per genome, which serves as the input to the consolidation stage described below.

## Motif Discovery and Consolidation

### STREME Motif Discovery
Motif discovery was performed using STREME (STochastic Regular Expression Motif Elicitation) version 5.4.1 from the MEME Suite (Bailey et al., 2015; Bailey, 2021). STREME was run independently on the promoter sequences of each genome using default parameters with the following settings: minimum motif width of 6 bp, maximum motif width of 15 bp, and a maximum of 20 motifs per analysis. The algorithm was configured to identify both palindromic and non-palindromic motifs using a first-order Markov model for background sequence generation.

### Motif Scanning with FIMO
To localise occurrences of the de novo motifs discovered by STREME, the FIMO program (Grant et al., 2011) from the MEME Suite was used to scan the same set of (masked) promoter sequences from which each genome's motifs were derived. This STREME → FIMO chain follows MEME-Suite recommendations: FIMO converts each motif into a log-odds position-specific scoring matrix and reports all sequence positions whose match score is significant under the configured p- or q-value threshold (defaults: p ≤ 1 × 10⁻⁴; alternatively a q-value cutoff of ≤ 0.05 is applied). Because the candidate promoters are ~1 kb in length, the expected per-promoter false-positive rate at the default threshold is low enough to keep biologically significant matches recoverable. A sequence-derived first-order Markov background model — the same model built for STREME — is supplied to FIMO via `--bgfile`, with fallback to the background embedded in the STREME motif file when an external background is not available. Scans can be run inline as part of the per-genome preparation pipeline or as a separate step, with one scan per genome dispatched in parallel.

### Motif Annotation with TOMTOM
Putative transcription-factor assignments for the de novo motifs were obtained with TOMTOM (Gupta et al., 2007) from the MEME Suite. Each genome's STREME motifs were compared against a reference plant TF position-weight-matrix database (e.g. PlantTFDB, JASPAR plants, or CIS-BP, all in MEME format), using the Pearson column-similarity metric and a default q-value threshold of 0.1 (configurable, with an E-value mode also available). TOMTOM was run per genome in parallel, with one `tomtom_<genome>/` output directory per genome containing the ranked match table and per-motif alignment plots. This step complements the FIMO occurrence analysis (which locates the motifs in the input sequences) by providing an external interpretation of what each de novo motif likely represents.

### Motif Consolidation Across Genomes
To address the challenge of comparing motifs discovered independently across multiple genomes, we implemented a comprehensive motif consolidation pipeline that accounts for the biological variability in motif instances while maintaining statistical rigor.

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
Gene expression data were processed to support both long format (Gene | Genome | Expression; the legacy `Line` column is accepted as an alias) and wide format (gene | LRTadd | one column per genome) inputs. For relative analysis, expression values represent log2 fold-changes relative to the designated reference genome (default: IM767, configurable via `--reference-genome`, with `--reference-line` retained as an alias), ensuring that the reference genome has expression values of zero across all genes.

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
The absolute analysis framework models expression as a direct function of motif features within each genome:

```
Expression_ij = β₀ + Σₖ(βₖ × Motif_Feature_ijk) + εᵢⱼ
```

Where:
- i = gene index
- j = genome index  
- k = motif feature index
- ε = residual error term

This approach is based on classical cis-regulatory theory where transcription factor binding sites additively contribute to gene expression (Bintu et al., 2005; Segal et al., 2008).

#### Relative Analysis (Novel Approach)
The relative analysis framework directly models expression differences as a function of regulatory feature differences from a reference genome:

```
ΔExpression_ij = β₀ + Σₖ(βₖ × ΔMotif_Feature_ijk) + εᵢⱼ
```

Where:
- ΔExpression_ij = Expression_ij - Expression_i,reference
- ΔMotif_Feature_ijk = Motif_Feature_ijk - Motif_Feature_i,reference,k

**Theoretical Justification**: This approach is grounded in evolutionary biology and comparative genomics, where regulatory evolution is understood through changes relative to ancestral states (Wittkopp & Kalay, 2012). By using a common reference genome (e.g. an isogenic reference line in single-species designs), we control for:
- Genetic background effects
- Trans-acting factors common across genomes
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
- **RepeatMasker** or **dust** (optional): Repeat/low-complexity masking of promoters

Promoter extraction from genome assemblies is implemented in pure Python (`cli_tools/genome_prep.py`) and requires no external dependencies. Masking, background modelling, and STREME are invoked through the same module, which detects the relevant external tool at run time.

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
1. **Independence**: Gene expression values are independent after accounting for genome
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