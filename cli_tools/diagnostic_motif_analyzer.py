#!/usr/bin/env python3
"""
Diagnostic Motif-Expression Analyzer

A comprehensive diagnostic tool to understand why motif-expression analysis
is producing poor R² values. This tool provides detailed statistics, 
visualizations, and recommendations for improving the analysis.

Author: Computational Biology Pipeline
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import GroupKFold, cross_val_score
from sklearn.linear_model import LinearRegression, Ridge, Lasso
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import r2_score, mean_squared_error
import warnings
import os
import argparse

warnings.filterwarnings('ignore')

def load_and_examine_data(motif_file, expr_file):
    """Load and provide detailed examination of input data"""
    print("="*70)
    print("DATA LOADING AND EXAMINATION")
    print("="*70)
    
    # Load motif data
    print(f"Loading motif data from: {motif_file}")
    motif_df = pd.read_csv(motif_file, sep='\t')
    print(f"Motif data shape: {motif_df.shape}")
    print(f"Columns: {list(motif_df.columns)}")
    print(f"Unique genes: {motif_df['gene_id'].nunique()}")
    print(f"Unique lines: {motif_df['line'].nunique()}")
    print(f"Unique motifs: {motif_df['consolidated_motif_id'].nunique()}")
    
    # Load expression data
    print(f"\nLoading expression data from: {expr_file}")
    expr_df = pd.read_csv(expr_file, sep='\t')
    print(f"Expression data shape: {expr_df.shape}")
    print(f"Columns: {list(expr_df.columns)}")
    
    # Check for wide format
    line_cols = [col for col in expr_df.columns if col != expr_df.columns[0]]
    print(f"Expression lines: {line_cols}")
    
    # Expression statistics
    expr_values = expr_df[line_cols].values.flatten()
    expr_values = expr_values[~np.isnan(expr_values)]
    
    print(f"\nExpression value statistics:")
    print(f"  Range: {expr_values.min():.3f} to {expr_values.max():.3f}")
    print(f"  Mean: {expr_values.mean():.3f}")
    print(f"  Std: {expr_values.std():.3f}")
    print(f"  Zero values: {(expr_values == 0).sum()} ({(expr_values == 0).mean()*100:.1f}%)")
    
    return motif_df, expr_df, line_cols

def analyze_motif_patterns(motif_df):
    """Analyze motif occurrence patterns"""
    print("\n" + "="*70)
    print("MOTIF PATTERN ANALYSIS")
    print("="*70)
    
    # Motif frequency per gene
    motif_per_gene = motif_df.groupby('gene_id')['consolidated_motif_id'].nunique()
    print(f"Motifs per gene - Mean: {motif_per_gene.mean():.1f}, "
          f"Median: {motif_per_gene.median():.1f}, "
          f"Range: {motif_per_gene.min()}-{motif_per_gene.max()}")
    
    # Motif frequency across lines
    motif_line_freq = motif_df.groupby(['consolidated_motif_id', 'line']).size().unstack(fill_value=0)
    motif_variance = motif_line_freq.var(axis=1)
    
    print(f"\nMotif line frequency variance:")
    print(f"  High variance motifs (top 10):")
    for motif, var in motif_variance.nlargest(10).items():
        print(f"    {motif}: {var:.3f}")
    
    print(f"  Low variance motifs (bottom 10):")
    for motif, var in motif_variance.nsmallest(10).items():
        print(f"    {motif}: {var:.3f}")
    
    return motif_variance

def create_feature_matrix(motif_df, expr_df, line_cols, reference_line='IM767'):
    """Create feature matrix with detailed diagnostics"""
    print("\n" + "="*70)
    print("FEATURE MATRIX CREATION")
    print("="*70)
    
    # Build motif features per gene-line
    motif_features = {}
    for _, row in motif_df.iterrows():
        gene = row['gene_id']
        line = row['line']
        motif = row['consolidated_motif_id']
        
        if gene not in motif_features:
            motif_features[gene] = {}
        if line not in motif_features[gene]:
            motif_features[gene][line] = {}
        if motif not in motif_features[gene][line]:
            motif_features[gene][line][motif] = {
                'present': 0, 'count': 0, 'positions': []
            }
        
        motif_features[gene][line][motif]['present'] = 1
        motif_features[gene][line][motif]['count'] += 1
        if 'start_pos' in row:
            motif_features[gene][line][motif]['positions'].append(row['start_pos'])
    
    print(f"Built motif features for {len(motif_features)} genes")
    
    # Get all motifs and comparison lines
    all_motifs = set()
    for gene_data in motif_features.values():
        for line_data in gene_data.values():
            all_motifs.update(line_data.keys())
    
    comparison_lines = [col for col in line_cols if col != reference_line]
    print(f"Total motifs: {len(all_motifs)}")
    print(f"Comparison lines: {len(comparison_lines)}")
    print(f"Reference line: {reference_line}")
    
    # Create relative features
    features_list = []
    
    for gene in expr_df[expr_df.columns[0]].values:
        gene_motifs = motif_features.get(gene, {})
        ref_motifs = gene_motifs.get(reference_line, {})
        
        for comp_line in comparison_lines:
            comp_motifs = gene_motifs.get(comp_line, {})
            
            # Get expression difference
            expr_diff = expr_df[expr_df[expr_df.columns[0]] == gene][comp_line].iloc[0]
            
            row_features = {
                'gene': gene,
                'comparison_line': comp_line,
                'expression_difference': expr_diff
            }
            
            # Add motif features
            for motif in all_motifs:
                ref_present = ref_motifs.get(motif, {}).get('present', 0)
                comp_present = comp_motifs.get(motif, {}).get('present', 0)
                
                ref_count = ref_motifs.get(motif, {}).get('count', 0)
                comp_count = comp_motifs.get(motif, {}).get('count', 0)
                
                row_features[f"{motif}_presence_diff"] = comp_present - ref_present
                row_features[f"{motif}_count_diff"] = comp_count - ref_count
            
            features_list.append(row_features)
    
    features_df = pd.DataFrame(features_list)
    print(f"Feature matrix shape: {features_df.shape}")
    
    return features_df

def diagnose_features(features_df):
    """Diagnose potential issues with features"""
    print("\n" + "="*70)
    print("FEATURE DIAGNOSTICS")
    print("="*70)
    
    feature_cols = [col for col in features_df.columns 
                   if col not in ['gene', 'comparison_line', 'expression_difference']]
    
    X = features_df[feature_cols]
    y = features_df['expression_difference']
    
    print(f"Feature matrix: {X.shape}")
    print(f"Expression differences: {y.shape}")
    
    # Check for constant features
    constant_features = X.columns[X.var() == 0]
    print(f"Constant features: {len(constant_features)}")
    
    # Check feature variance
    feature_vars = X.var().sort_values(ascending=False)
    print(f"\nTop 10 highest variance features:")
    for feat, var in feature_vars.head(10).items():
        print(f"  {feat}: {var:.6f}")
    
    print(f"\nBottom 10 lowest variance features:")
    for feat, var in feature_vars.tail(10).items():
        print(f"  {feat}: {var:.6f}")
    
    # Check for non-zero features
    non_zero_counts = (X != 0).sum()
    print(f"\nFeatures with most non-zero values:")
    for feat, count in non_zero_counts.sort_values(ascending=False).head(10).items():
        print(f"  {feat}: {count} ({count/len(X)*100:.1f}%)")
    
    # Correlation with target
    print(f"\nExpression difference statistics:")
    print(f"  Range: {y.min():.3f} to {y.max():.3f}")
    print(f"  Mean: {y.mean():.3f}")
    print(f"  Std: {y.std():.3f}")
    print(f"  Zeros: {(y == 0).sum()} ({(y == 0).mean()*100:.1f}%)")
    
    # Feature-target correlations
    correlations = []
    for col in feature_cols:
        if X[col].var() > 0:  # Skip constant features
            corr = np.corrcoef(X[col], y)[0, 1]
            if not np.isnan(corr):
                correlations.append((col, abs(corr)))
    
    correlations.sort(key=lambda x: x[1], reverse=True)
    
    print(f"\nTop 10 features by absolute correlation with expression:")
    for feat, corr in correlations[:10]:
        print(f"  {feat}: {corr:.4f}")
    
    return X, y, feature_vars, correlations

def test_simple_models(X, y, features_df):
    """Test simple models with diagnostics"""
    print("\n" + "="*70)
    print("SIMPLE MODEL TESTING")
    print("="*70)
    
    # Remove constant features
    constant_features = X.columns[X.var() == 0]
    if len(constant_features) > 0:
        print(f"Removing {len(constant_features)} constant features")
        X = X.drop(columns=constant_features)
    
    # Test different preprocessing approaches
    approaches = {
        'Raw Features': X,
        'Standardized': pd.DataFrame(StandardScaler().fit_transform(X), 
                                   columns=X.columns, index=X.index),
    }
    
    # Simple models
    models = {
        'Linear Regression': LinearRegression(),
        'Ridge (alpha=1.0)': Ridge(alpha=1.0),
        'Ridge (alpha=10.0)': Ridge(alpha=10.0),
    }
    
    results = []
    
    for approach_name, X_processed in approaches.items():
        print(f"\n{approach_name}:")
        
        for model_name, model in models.items():
            # Different CV strategies
            cv_strategies = {
                'Standard 5-Fold': 5,
                'Gene Groups': GroupKFold(n_splits=5)
            }
            
            for cv_name, cv in cv_strategies.items():
                try:
                    if cv_name == 'Gene Groups':
                        groups = features_df['gene']
                        cv_scores = cross_val_score(model, X_processed, y, 
                                                  cv=cv, groups=groups, scoring='r2')
                    else:
                        cv_scores = cross_val_score(model, X_processed, y, 
                                                  cv=cv, scoring='r2')
                    
                    # Full model fit
                    model.fit(X_processed, y)
                    y_pred = model.predict(X_processed)
                    full_r2 = r2_score(y, y_pred)
                    
                    print(f"  {model_name:20} ({cv_name:15}): "
                          f"CV R² = {cv_scores.mean():.4f} ± {cv_scores.std():.4f}, "
                          f"Full R² = {full_r2:.4f}")
                    
                    results.append({
                        'approach': approach_name,
                        'model': model_name,
                        'cv_strategy': cv_name,
                        'cv_r2_mean': cv_scores.mean(),
                        'cv_r2_std': cv_scores.std(),
                        'full_r2': full_r2
                    })
                    
                except Exception as e:
                    print(f"  {model_name:20} ({cv_name:15}): Error - {e}")
    
    return pd.DataFrame(results)

def generate_recommendations(results_df, correlations, feature_vars):
    """Generate recommendations based on diagnostics"""
    print("\n" + "="*70)
    print("RECOMMENDATIONS")
    print("="*70)
    
    best_result = results_df.loc[results_df['cv_r2_mean'].idxmax()]
    print(f"Best performing setup:")
    print(f"  {best_result['approach']} + {best_result['model']} + {best_result['cv_strategy']}")
    print(f"  CV R² = {best_result['cv_r2_mean']:.4f}")
    
    print("\nDiagnostic insights:")
    
    # Check if any correlations exist
    max_corr = max(correlations, key=lambda x: x[1])[1] if correlations else 0
    print(f"  Highest feature-expression correlation: {max_corr:.4f}")
    
    if max_corr < 0.01:
        print("  ⚠️  Very weak correlations - features may not be predictive")
    
    # Check feature variance
    high_var_features = (feature_vars > 0.1).sum()
    print(f"  Features with variance > 0.1: {high_var_features}")
    
    if high_var_features < 10:
        print("  ⚠️  Few high-variance features - consider different feature engineering")
    
    print("\nRecommendations:")
    
    if max_corr < 0.05:
        print("1. 🔧 Consider different feature engineering approaches:")
        print("   - Try absolute analysis instead of relative")
        print("   - Focus on specific motif types or gene subsets")
        print("   - Include interaction terms or non-linear features")
    
    if high_var_features < 20:
        print("2. 📊 Improve feature selection:")
        print("   - Try different motif selection criteria")
        print("   - Consider motif families or functional groups")
        print("   - Use domain knowledge to filter motifs")
    
    if best_result['cv_r2_mean'] < 0.01:
        print("3. 🔬 Consider alternative analysis approaches:")
        print("   - Gene-specific analysis for individual genes")
        print("   - Pathway-based analysis")
        print("   - Non-linear modeling approaches")
    
    print("4. 🧪 Experimental validation:")
    print("   - Check if expression differences are real (not technical noise)")
    print("   - Validate motif predictions with known TF binding sites")
    print("   - Consider chromatin accessibility data")

def main():
    parser = argparse.ArgumentParser(description="Diagnostic motif-expression analysis")
    parser.add_argument('motif_file', help='Consolidated motif sites file')
    parser.add_argument('expr_file', help='Expression data file')
    parser.add_argument('--reference-line', default='IM767', 
                       help='Reference line for relative analysis')
    parser.add_argument('--output', default='diagnostic_results', 
                       help='Output directory')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    # Run diagnostics
    motif_df, expr_df, line_cols = load_and_examine_data(args.motif_file, args.expr_file)
    motif_variance = analyze_motif_patterns(motif_df)
    features_df = create_feature_matrix(motif_df, expr_df, line_cols, args.reference_line)
    X, y, feature_vars, correlations = diagnose_features(features_df)
    results_df = test_simple_models(X, y, features_df)
    generate_recommendations(results_df, correlations, feature_vars)
    
    # Save results
    results_file = os.path.join(args.output, 'diagnostic_results.csv')
    results_df.to_csv(results_file, index=False)
    print(f"\n💾 Diagnostic results saved to: {results_file}")
    
    # Save feature correlations
    corr_file = os.path.join(args.output, 'feature_correlations.csv')
    pd.DataFrame(correlations, columns=['Feature', 'Abs_Correlation']).to_csv(corr_file, index=False)
    print(f"💾 Feature correlations saved to: {corr_file}")
    
    print(f"\n🎉 Diagnostic analysis complete!")

if __name__ == "__main__":
    main()