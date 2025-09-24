#!/usr/bin/env python3
"""
Relative Motif Analyzer - Compare regulatory landscapes against reference baseline

This tool performs differential regulatory analysis by comparing each line's
motif profile against a reference line (e.g., IM767) to identify regulatory
changes that drive expression differences.

Key Features:
- Uses IM767 (or specified line) as baseline reference (expression = 0)
- Generates relative regulatory features (LineX vs Reference)
- Models how regulatory differences drive expression changes
- Focuses on gained/lost motifs, position shifts, sequence variants
"""

import os
import sys
import pandas as pd
import numpy as np
import argparse
from pathlib import Path
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.linear_model import LinearRegression, LassoCV, RidgeCV
from sklearn.model_selection import cross_val_score, GroupKFold
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.preprocessing import StandardScaler
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

def load_consolidated_motif_data(motif_file):
    """Load consolidated STREME sites data"""
    print(f"Loading consolidated motif data from: {motif_file}")
    
    df = pd.read_csv(motif_file, sep='\t')
    
    # Ensure required columns exist
    required_cols = ['gene_id', 'line', 'consolidated_motif_id', 'motif_consensus', 'sequence', 'start_pos']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"ERROR: Missing required columns: {missing_cols}")
        return None
    
    print(f"Loaded {len(df)} motif sites")
    print(f"Genes: {df['gene_id'].nunique()}")
    print(f"Lines: {df['line'].nunique()}")
    print(f"Motifs: {df['consolidated_motif_id'].nunique()}")
    
    return df

def load_wide_expression_data(expression_file, reference_line='IM767'):
    """
    Load wide-format expression data where genes are rows and lines are columns
    
    Expected format:
    gene    LRTadd  IM62    IM155   IM444   IM502   IM541   IM664   IM909   IM1034  IM1192  IM767
    MiIM7v11000007m.g  47.272  0.332   0.078   0.005   0.074   0.294   0.443   0.458   0.1     0.11    0.0
    """
    print(f"Loading wide-format expression data from: {expression_file}")
    
    df = pd.read_csv(expression_file, sep='\t')
    
    # Get gene column (first column)
    gene_col = df.columns[0]
    
    # Get line columns (exclude gene and LRTadd if present)
    line_cols = [col for col in df.columns if col not in [gene_col, 'LRTadd']]
    
    # Verify reference line exists
    if reference_line not in line_cols:
        print(f"ERROR: Reference line {reference_line} not found in columns: {line_cols}")
        return None, None
    
    # Check that reference line has all zeros (or close to zero)
    ref_values = df[reference_line].values
    if not np.allclose(ref_values, 0.0, atol=1e-6):
        print(f"WARNING: Reference line {reference_line} values are not all zero")
        print(f"Range: {ref_values.min()} to {ref_values.max()}")
    
    print(f"Loaded expression data for {len(df)} genes across {len(line_cols)} lines")
    print(f"Reference line: {reference_line}")
    print(f"Comparison lines: {[col for col in line_cols if col != reference_line]}")
    
    return df, line_cols

def calculate_motif_features_per_line(motif_df):
    """
    Calculate motif features for each gene×line combination
    Returns: dict[gene][line][motif] = features_dict
    """
    print("Calculating motif features per gene-line combination...")
    
    features = defaultdict(lambda: defaultdict(lambda: defaultdict(dict)))
    
    # Group by gene, line, motif
    for (gene, line, motif), group in motif_df.groupby(['gene_id', 'line', 'consolidated_motif_id']):
        positions = group['start_pos'].values
        sequences = group['sequence'].values
        
        # Basic presence and count
        features[gene][line][motif]['present'] = 1
        features[gene][line][motif]['count'] = len(group)
        
        # Position features
        features[gene][line][motif]['position_mean'] = np.mean(positions)
        features[gene][line][motif]['position_std'] = np.std(positions) if len(positions) > 1 else 0
        features[gene][line][motif]['position_min'] = np.min(positions)
        features[gene][line][motif]['position_max'] = np.max(positions)
        
        # Position density (assuming 1000bp promoter)
        proximal_count = np.sum((positions >= 700) & (positions <= 1000))  # Close to TSS
        distal_count = np.sum((positions >= 0) & (positions < 300))        # Far from TSS
        
        features[gene][line][motif]['proximal_density'] = proximal_count / 3.0  # per 100bp
        features[gene][line][motif]['distal_density'] = distal_count / 3.0
        
        # Sequence variation features
        unique_seqs = len(set(sequences))
        features[gene][line][motif]['sequence_diversity'] = unique_seqs / len(sequences)
        
        # GC content
        gc_contents = []
        for seq in sequences:
            gc_count = seq.upper().count('G') + seq.upper().count('C')
            gc_content = gc_count / len(seq) if seq else 0
            gc_contents.append(gc_content)
        
        features[gene][line][motif]['gc_content_mean'] = np.mean(gc_contents)
        features[gene][line][motif]['gc_content_std'] = np.std(gc_contents) if len(gc_contents) > 1 else 0
    
    return features

def generate_relative_features(motif_features, expr_df, line_cols, reference_line='IM767', 
                             top_motifs=None, min_sites=10, selection_method='frequency'):
    """
    Generate relative regulatory features comparing each line to the reference line
    
    selection_method options:
    - 'frequency': Most common motifs (default, conservative)
    - 'variance': Motifs with highest presence variance across lines (differential)
    - 'expression_corr': Motifs with highest correlation to expression differences
    """
    print(f"Generating relative features with {reference_line} as baseline...")
    
    # Get all motifs and calculate statistics
    all_motifs = set()
    motif_counts = defaultdict(int)
    motif_line_presence = defaultdict(lambda: defaultdict(int))
    
    for gene, gene_features in motif_features.items():
        for line, line_features in gene_features.items():
            for motif in line_features.keys():
                all_motifs.add(motif)
                motif_counts[motif] += 1
                motif_line_presence[motif][line] += 1
    
    # Filter motifs by minimum frequency
    if min_sites:
        all_motifs = {motif for motif in all_motifs if motif_counts[motif] >= min_sites}
        print(f"Filtered to {len(all_motifs)} motifs with >= {min_sites} sites")
    
    if top_motifs and top_motifs < len(all_motifs):
        if selection_method == 'frequency':
            # Select top motifs by frequency (original method)
            sorted_motifs = sorted(motif_counts.items(), key=lambda x: x[1], reverse=True)
            all_motifs = {motif for motif, count in sorted_motifs[:top_motifs]}
            print(f"Selected top {top_motifs} most frequent motifs")
            
        elif selection_method == 'variance':
            # Select motifs with highest presence variance across lines
            motif_variances = {}
            comparison_lines = [col for col in line_cols if col != reference_line]
            
            for motif in all_motifs:
                # Calculate presence frequency per line
                line_frequencies = []
                for line in comparison_lines:
                    total_genes = len([g for g in motif_features.keys() 
                                     if line in motif_features[g]])
                    present_genes = len([g for g in motif_features.keys() 
                                       if line in motif_features[g] and motif in motif_features[g][line]])
                    freq = present_genes / total_genes if total_genes > 0 else 0
                    line_frequencies.append(freq)
                
                # Calculate variance in presence across lines
                motif_variances[motif] = np.var(line_frequencies) if len(line_frequencies) > 1 else 0
            
            # Select top variable motifs
            sorted_motifs = sorted(motif_variances.items(), key=lambda x: x[1], reverse=True)
            all_motifs = {motif for motif, var in sorted_motifs[:top_motifs]}
            print(f"Selected top {top_motifs} most variable motifs (presence variance)")
            
        elif selection_method == 'expression_corr':
            # Select motifs most correlated with expression differences
            # This is more complex - implement a simplified version
            print("Warning: expression_corr selection not fully implemented, using variance method")
            # Fall back to variance method
            motif_variances = {}
            comparison_lines = [col for col in line_cols if col != reference_line]
            
            for motif in all_motifs:
                line_frequencies = []
                for line in comparison_lines:
                    total_genes = len([g for g in motif_features.keys() 
                                     if line in motif_features[g]])
                    present_genes = len([g for g in motif_features.keys() 
                                       if line in motif_features[g] and motif in motif_features[g][line]])
                    freq = present_genes / total_genes if total_genes > 0 else 0
                    line_frequencies.append(freq)
                
                motif_variances[motif] = np.var(line_frequencies) if len(line_frequencies) > 1 else 0
            
            sorted_motifs = sorted(motif_variances.items(), key=lambda x: x[1], reverse=True)
            all_motifs = {motif for motif, var in sorted_motifs[:top_motifs]}
            print(f"Selected top {top_motifs} most variable motifs")
    
    all_motifs = sorted(all_motifs)
    comparison_lines = [col for col in line_cols if col != reference_line]
    
    # Generate relative features for each gene×comparison_line
    relative_data = []
    
    for gene in expr_df[expr_df.columns[0]].values:
        gene_motif_features = motif_features.get(gene, {})
        ref_features = gene_motif_features.get(reference_line, {})
        
        for comp_line in comparison_lines:
            comp_features = gene_motif_features.get(comp_line, {})
            
            # Get expression difference (already computed in the data)
            expr_diff = expr_df[expr_df[expr_df.columns[0]] == gene][comp_line].iloc[0]
            
            # Generate relative features for each motif
            row_features = {
                'gene': gene,
                'comparison_line': comp_line,
                'expression_difference': expr_diff
            }
            
            for motif in all_motifs:
                ref_motif = ref_features.get(motif, {})
                comp_motif = comp_features.get(motif, {})
                
                # Presence difference (gained/lost motif)
                ref_present = ref_motif.get('present', 0)
                comp_present = comp_motif.get('present', 0)
                row_features[f"{motif}_presence_diff"] = comp_present - ref_present
                
                # Count difference
                ref_count = ref_motif.get('count', 0)
                comp_count = comp_motif.get('count', 0)
                row_features[f"{motif}_count_diff"] = comp_count - ref_count
                
                # Position differences (only if both lines have the motif)
                if ref_present and comp_present:
                    ref_pos_mean = ref_motif.get('position_mean', 0)
                    comp_pos_mean = comp_motif.get('position_mean', 0)
                    row_features[f"{motif}_position_mean_diff"] = comp_pos_mean - ref_pos_mean
                    
                    # Density differences
                    ref_prox = ref_motif.get('proximal_density', 0)
                    comp_prox = comp_motif.get('proximal_density', 0)
                    row_features[f"{motif}_proximal_density_diff"] = comp_prox - ref_prox
                    
                    # Sequence diversity differences
                    ref_div = ref_motif.get('sequence_diversity', 0)
                    comp_div = comp_motif.get('sequence_diversity', 0)
                    row_features[f"{motif}_sequence_diversity_diff"] = comp_div - ref_div
                    
                    # GC content differences
                    ref_gc = ref_motif.get('gc_content_mean', 0)
                    comp_gc = comp_motif.get('gc_content_mean', 0)
                    row_features[f"{motif}_gc_content_diff"] = comp_gc - ref_gc
                    
                else:
                    # Set to 0 if motif not present in both lines
                    row_features[f"{motif}_position_mean_diff"] = 0
                    row_features[f"{motif}_proximal_density_diff"] = 0
                    row_features[f"{motif}_sequence_diversity_diff"] = 0
                    row_features[f"{motif}_gc_content_diff"] = 0
            
            relative_data.append(row_features)
    
    relative_df = pd.DataFrame(relative_data)
    
    print(f"Generated relative features: {len(relative_df)} gene×line comparisons")
    print(f"Feature matrix shape: {relative_df.shape}")
    
    return relative_df

def run_relative_analysis(relative_df):
    """
    Run regression analysis on relative features
    """
    print("\nRunning relative regulatory analysis...")
    
    # Prepare feature matrix
    feature_cols = [col for col in relative_df.columns 
                   if col not in ['gene', 'comparison_line', 'expression_difference']]
    
    X = relative_df[feature_cols]
    y = relative_df['expression_difference']
    
    # Remove constant features
    constant_features = X.columns[X.var() == 0]
    if len(constant_features) > 0:
        print(f"Removing {len(constant_features)} constant features")
        X = X.drop(columns=constant_features)
    
    print(f"Final feature matrix: {X.shape[0]} samples × {X.shape[1]} features")
    print(f"Expression difference range: {y.min():.3f} to {y.max():.3f}")
    
    # Use GroupKFold to ensure same gene doesn't appear in train and test
    cv = GroupKFold(n_splits=5)
    groups = relative_df['gene']
    
    # Test multiple models
    models = {
        'Linear Regression': LinearRegression(),
        'Ridge Regression': RidgeCV(alphas=[0.1, 1.0, 10.0, 100.0]),
        'Lasso Regression': LassoCV(cv=3, max_iter=10000),
        'Random Forest': RandomForestRegressor(n_estimators=100, random_state=42),
        'Gradient Boosting': GradientBoostingRegressor(n_estimators=100, random_state=42)
    }
    
    results = {}
    
    for name, model in models.items():
        try:
            # Cross-validation scores
            cv_scores = cross_val_score(model, X, y, cv=cv, groups=groups, 
                                      scoring='r2', n_jobs=-1)
            
            # Fit full model for feature importance
            model.fit(X, y)
            y_pred = model.predict(X)
            
            results[name] = {
                'cv_r2_mean': cv_scores.mean(),
                'cv_r2_std': cv_scores.std(),
                'full_r2': r2_score(y, y_pred),
                'rmse': np.sqrt(mean_squared_error(y, y_pred)),
                'model': model
            }
            
            print(f"{name:20s} - CV R²: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
            
        except Exception as e:
            print(f"Error with {name}: {e}")
            continue
    
    return results, X, y

def analyze_relative_feature_importance(results, X, top_n=20):
    """Analyze which regulatory differences are most important"""
    print(f"\nAnalyzing top {top_n} most important regulatory differences...")
    
    importance_data = []
    
    for model_name, result in results.items():
        model = result['model']
        
        if hasattr(model, 'feature_importances_'):
            # Tree-based models
            importances = model.feature_importances_
        elif hasattr(model, 'coef_'):
            # Linear models
            importances = np.abs(model.coef_)
        else:
            continue
        
        # Get top features
        feature_names = X.columns
        importance_df = pd.DataFrame({
            'Feature': feature_names,
            'Importance': importances,
            'Model': model_name
        }).sort_values('Importance', ascending=False)
        
        importance_data.append(importance_df.head(top_n))
    
    if importance_data:
        all_importance = pd.concat(importance_data, ignore_index=True)
        return all_importance
    else:
        return pd.DataFrame()

def create_relative_visualizations(results, importance_df, relative_df, output_dir):
    """Create visualizations for relative analysis"""
    print("Creating relative analysis visualizations...")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Model performance comparison
    plt.figure(figsize=(10, 6))
    model_names = list(results.keys())
    cv_scores = [results[name]['cv_r2_mean'] for name in model_names]
    cv_errors = [results[name]['cv_r2_std'] for name in model_names]
    
    plt.bar(model_names, cv_scores, yerr=cv_errors, capsize=5)
    plt.ylabel('Cross-Validation R²')
    plt.title('Relative Regulatory Analysis - Model Performance')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/relative_model_performance.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Feature importance heatmap
    if not importance_df.empty:
        plt.figure(figsize=(14, 10))
        pivot_importance = importance_df.pivot(index='Feature', columns='Model', values='Importance')
        sns.heatmap(pivot_importance.head(20), annot=True, fmt='.3f', cmap='RdYlBu_r')
        plt.title('Top 20 Regulatory Differences - Feature Importance')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/relative_feature_importance_heatmap.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    # 3. Expression difference distribution by line
    plt.figure(figsize=(12, 6))
    sns.boxplot(data=relative_df, x='comparison_line', y='expression_difference')
    plt.ylabel('Expression Difference from IM767')
    plt.xlabel('Genetic Line')
    plt.title('Expression Difference Distribution by Line (vs IM767)')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/expression_differences_by_line.png", dpi=300, bbox_inches='tight')
    plt.close()

def generate_relative_report(results, importance_df, relative_df, output_file, reference_line='IM767'):
    """Generate comprehensive relative analysis report"""
    print(f"Generating relative analysis report: {output_file}")
    
    with open(output_file, 'w') as f:
        f.write("# Relative Motif-Expression Analysis Report\n\n")
        f.write(f"**Analysis Type**: Differential regulatory analysis vs {reference_line} baseline\n\n")
        
        # Data summary
        f.write("## Data Summary\n")
        f.write(f"- Total gene×line comparisons: {len(relative_df)}\n")
        f.write(f"- Unique genes: {relative_df['gene'].nunique()}\n")
        f.write(f"- Reference line: {reference_line}\n")
        f.write(f"- Comparison lines: {', '.join(sorted(relative_df['comparison_line'].unique()))}\n")
        f.write(f"- Expression difference range: {relative_df['expression_difference'].min():.3f} to {relative_df['expression_difference'].max():.3f}\n\n")
        
        # Model performance
        f.write("## Model Performance\n")
        f.write("| Model | CV R² | CV Std | Full R² | RMSE |\n")
        f.write("|-------|-------|--------|---------|------|\n")
        for name, result in results.items():
            f.write(f"| {name} | {result['cv_r2_mean']:.3f} | {result['cv_r2_std']:.3f} | "
                   f"{result['full_r2']:.3f} | {result['rmse']:.3f} |\n")
        f.write("\n")
        
        # Best model
        best_model = max(results.keys(), key=lambda k: results[k]['cv_r2_mean'])
        f.write(f"**Best performing model:** {best_model} (CV R² = {results[best_model]['cv_r2_mean']:.3f})\n\n")
        
        # Top regulatory differences
        if not importance_df.empty:
            f.write("## Top 10 Most Important Regulatory Differences\n")
            top_features = importance_df.groupby('Feature')['Importance'].mean().sort_values(ascending=False).head(10)
            for feature, importance in top_features.items():
                # Parse feature name
                if '_presence_diff' in feature:
                    motif = feature.replace('_presence_diff', '')
                    feature_type = "Presence Change"
                elif '_count_diff' in feature:
                    motif = feature.replace('_count_diff', '')
                    feature_type = "Count Change"
                elif '_position_mean_diff' in feature:
                    motif = feature.replace('_position_mean_diff', '')
                    feature_type = "Position Shift"
                elif '_proximal_density_diff' in feature:
                    motif = feature.replace('_proximal_density_diff', '')
                    feature_type = "Proximal Density Change"
                else:
                    motif = feature
                    feature_type = "Other"
                
                f.write(f"- **{motif}** ({feature_type}): {importance:.3f}\n")
            f.write("\n")
        
        # Interpretation
        f.write("## Interpretation\n")
        r2_best = results[best_model]['cv_r2_mean']
        if r2_best > 0.4:
            f.write("🎉 **Strong predictive power**: Regulatory differences vs baseline explain substantial expression variation.\n")
        elif r2_best > 0.2:
            f.write("📈 **Moderate predictive power**: Regulatory differences contribute to expression changes.\n")
        else:
            f.write("⚠️ **Limited predictive power**: Consider other regulatory mechanisms or data quality.\n")
        
        f.write(f"\nThe analysis suggests that regulatory differences relative to {reference_line} "
               f"can explain {r2_best*100:.1f}% of the variance in expression differences.\n\n")
        
        f.write("## Key Insights\n")
        f.write("- **Gained/Lost Motifs**: Which motifs, when gained or lost relative to baseline, drive expression changes?\n")
        f.write("- **Position Shifts**: How do changes in motif position relative to TSS affect expression?\n")
        f.write("- **Copy Number Changes**: Does having more/fewer copies of a motif change expression?\n")
        f.write("- **Sequence Variation**: Do different variants of the same motif family have different effects?\n")

def main():
    parser = argparse.ArgumentParser(
        description='Relative motif analysis - compare regulatory landscapes vs reference baseline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic relative analysis vs IM767
  %(prog)s consolidated_sites.tsv expression_wide.tsv --output results/relative_analysis/
  
  # Use different reference line
  %(prog)s consolidated_sites.tsv expression_wide.tsv --reference-line IM502 --output results/
  
  # Focus on top motifs
  %(prog)s consolidated_sites.tsv expression_wide.tsv --top-motifs 50 --min-sites 20

Expected expression format (wide):
  gene    LRTadd  IM62    IM155   ...  IM767
  Gene1   47.272  0.332   0.078   ...  0.0
  Gene2   23.36   -0.05   -0.108  ...  0.0
        """
    )
    
    parser.add_argument('motif_sites_file', help='Consolidated STREME sites TSV file')
    parser.add_argument('expression_data', help='Wide-format expression data file')
    parser.add_argument('--output', '-o', default='relative_analysis_results', 
                       help='Output directory for results')
    parser.add_argument('--reference-line', '-r', default='IM767',
                       help='Reference line to use as baseline (default: IM767)')
    parser.add_argument('--top-motifs', '-t', type=int, default=100,
                       help='Number of top motifs to include (default: 100)')
    parser.add_argument('--min-sites', '-m', type=int, default=10,
                       help='Minimum sites required for motif inclusion (default: 10)')
    parser.add_argument('--selection-method', choices=['frequency', 'variance', 'expression_corr'], 
                       default='variance',
                       help='Method for selecting top motifs: frequency (most common), variance (most differential), expression_corr (correlated with expression)')
    
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    print("=" * 70)
    print("RELATIVE MOTIF-EXPRESSION ANALYSIS")
    print(f"Reference baseline: {args.reference_line}")
    print("=" * 70)
    
    try:
        # Load consolidated motif data
        motif_df = load_consolidated_motif_data(args.motif_sites_file)
        if motif_df is None:
            return 1
        
        # Load wide-format expression data
        expr_df, line_cols = load_wide_expression_data(args.expression_data, args.reference_line)
        if expr_df is None:
            return 1
        
        # Calculate motif features per line
        motif_features = calculate_motif_features_per_line(motif_df)
        
        # Generate relative features
        relative_df = generate_relative_features(
            motif_features, expr_df, line_cols, args.reference_line, 
            args.top_motifs, args.min_sites, args.selection_method
        )
        
        if len(relative_df) == 0:
            print("ERROR: No relative features generated.")
            return 1
        
        # Run relative analysis
        results, X, y = run_relative_analysis(relative_df)
        
        if not results:
            print("ERROR: No models ran successfully")
            return 1
        
        # Feature importance analysis
        importance_df = analyze_relative_feature_importance(results, X)
        
        # Create visualizations
        create_relative_visualizations(results, importance_df, relative_df, args.output)
        
        # Generate report
        report_file = os.path.join(args.output, 'relative_analysis_report.md')
        generate_relative_report(results, importance_df, relative_df, report_file, args.reference_line)
        
        # Save detailed results
        if not importance_df.empty:
            importance_df.to_csv(os.path.join(args.output, 'relative_feature_importance.tsv'), 
                                sep='\t', index=False)
        
        # Save relative features and predictions
        relative_df.to_csv(os.path.join(args.output, 'relative_features.tsv'), 
                          sep='\t', index=False)
        
        # Save predictions from best model
        best_model_name = max(results.keys(), key=lambda k: results[k]['cv_r2_mean'])
        best_model = results[best_model_name]['model']
        
        predictions_df = relative_df[['gene', 'comparison_line', 'expression_difference']].copy()
        predictions_df['predicted_expression_difference'] = best_model.predict(X)
        predictions_df['residual'] = predictions_df['expression_difference'] - predictions_df['predicted_expression_difference']
        predictions_df.to_csv(os.path.join(args.output, 'relative_predictions.tsv'), 
                             sep='\t', index=False)
        
        print(f"\n🎉 Relative analysis complete! Results saved to: {args.output}/")
        print(f"📊 Check relative_analysis_report.md for detailed findings")
        print(f"📈 Best model: {best_model_name} (R² = {results[best_model_name]['cv_r2_mean']:.3f})")
        print(f"🔍 This explains {results[best_model_name]['cv_r2_mean']*100:.1f}% of expression differences vs {args.reference_line}")
        
        return 0
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())