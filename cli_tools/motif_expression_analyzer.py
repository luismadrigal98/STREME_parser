#!/usr/bin/env python3
"""
Absolute Motif-Expression Analysis Tool

This tool performs integrative analysis of motif occurrences and gene expression data
to identify regulatory relationships and patterns. This provides ABSOLUTE analysis
(per-line analysis). For RELATIVE analysis (comparing to IM767 baseline), 
use relative_motif_analyzer.py instead.

Features:
- Load consolidated motif data from STREME analysis
- Integrate with gene expression data (long or wide format)
- Generate comprehensive regulatory feature matrices
- Perform correlation analysis between motifs and expression
- Model gene expression using motif features
- Generate publication-ready visualizations and reports

Usage:
    python motif_expression_analyzer.py <consolidated_motifs.tsv> <expression_data.tsv> [options]

Input formats:
    - Motif data: Output from streme_sites_consolidator.py
    - Expression data: 
      * Long format: Gene | Line | Expression (tab-separated)
      * Wide format: gene | LRTadd | IM62 | IM155 | ... | IM767 (tab-separated)

Output:
    - Regulatory feature matrix
    - Expression correlation analysis
    - Regression models and coefficients
    - Comprehensive visualizations
    - Summary statistics and reports

Author: Computational Biology Pipeline
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
import warnings
warnings.filterwarnings('ignore')

def load_motif_features(motif_file):
    """Load motif features from extract-features output"""
    print(f"Loading motif features from: {motif_file}")
    
    # Load the feature matrix
    df = pd.read_csv(motif_file, sep='\t')
    
    print(f"Found {len(df)} gene-line combinations")
    print(f"Found {len([c for c in df.columns if c.endswith('_present')])} motifs")
    
    return df

def load_expression_data(expression_file):
    """
    Load expression data. Supports two formats:
    1. Long format: Gene | Line | Expression
    2. Wide format: gene | LRTadd | IM62 | IM155 | ... | IM767
    """
    print(f"Loading expression data from: {expression_file}")
    
    expr_df = pd.read_csv(expression_file, sep='\t')
    
    # Detect format based on columns
    if len(expr_df.columns) > 3 and any(col.startswith('IM') for col in expr_df.columns):
        # Wide format - convert to long format
        print("Detected wide format - converting to long format")
        
        gene_col = expr_df.columns[0]
        line_cols = [col for col in expr_df.columns if col not in [gene_col, 'LRTadd'] and col.startswith('IM')]
        
        # Melt to long format
        long_df = expr_df.melt(
            id_vars=[gene_col], 
            value_vars=line_cols,
            var_name='Line', 
            value_name='Expression'
        )
        long_df = long_df.rename(columns={gene_col: 'Gene'})
        
        # Remove IM767 baseline (expression = 0)
        long_df = long_df[long_df['Line'] != 'IM767'].copy()
        
        print(f"Converted wide format to long format")
        print(f"Found expression data for {long_df['Gene'].nunique()} genes across {long_df['Line'].nunique()} lines")
        print(f"Lines found: {sorted(long_df['Line'].unique())}")
        
        return long_df
    
    else:
        # Long format (original)
        # Ensure proper column names
        expected_cols = ['Gene', 'Line', 'Expression']
        if list(expr_df.columns) != expected_cols:
            print(f"Warning: Expected columns {expected_cols}, got {list(expr_df.columns)}")
            print("Assuming first 3 columns are Gene, Line, Expression")
            expr_df.columns = expected_cols[:len(expr_df.columns)]
        
        # Remove IM767 entries (they should be 0 anyway)
        expr_df = expr_df[expr_df['Line'] != 'IM767'].copy()
        
        print(f"Found expression data for {expr_df['Gene'].nunique()} genes across {expr_df['Line'].nunique()} lines")
        print(f"Lines found: {sorted(expr_df['Line'].unique())}")
        
        return expr_df

def merge_motif_expression_data(motif_df, expr_df):
    """Merge motif features with expression data"""
    print("Merging motif features with expression data...")
    
    # Ensure consistent naming
    if 'gene_id' in motif_df.columns:
        motif_df = motif_df.rename(columns={'gene_id': 'Gene'})
    if 'line' in motif_df.columns:
        motif_df = motif_df.rename(columns={'line': 'Line'})
    
    # Merge on Gene and Line
    merged_df = pd.merge(motif_df, expr_df, on=['Gene', 'Line'], how='inner')
    
    print(f"Successfully merged data for {len(merged_df)} gene-line combinations")
    print(f"Genes with both motif and expression data: {merged_df['Gene'].nunique()}")
    
    return merged_df

def prepare_regression_data(merged_df, simple_mode=True):
    """Prepare data for regression analysis"""
    print("Preparing regression data...")
    
    # Identify motif columns
    if simple_mode:
        motif_cols = [col for col in merged_df.columns if col.endswith('_present')]
    else:
        motif_cols = [col for col in merged_df.columns 
                     if not col in ['Gene', 'Line', 'Expression'] and 
                        any(col.endswith(suffix) for suffix in ['_present', '_count', '_mean', '_std'])]
    
    print(f"Using {len(motif_cols)} motif features")
    
    # Prepare feature matrix (X) and target (y)
    X = merged_df[motif_cols].copy()
    y = merged_df['Expression'].copy()
    
    # Additional metadata for analysis
    metadata = merged_df[['Gene', 'Line']].copy()
    
    # Remove constant features (motifs present in all or no samples)
    constant_features = X.columns[X.var() == 0]
    if len(constant_features) > 0:
        print(f"Removing {len(constant_features)} constant features")
        X = X.drop(columns=constant_features)
    
    print(f"Final feature matrix: {X.shape[0]} samples × {X.shape[1]} features")
    print(f"Expression range: {y.min():.2f} to {y.max():.2f}")
    
    return X, y, metadata, motif_cols

def run_regression_analysis(X, y, metadata):
    """Run multiple regression models and compare performance"""
    print("\nRunning regression analysis...")
    
    # Use GroupKFold to ensure same gene doesn't appear in train and test
    cv = GroupKFold(n_splits=5)
    groups = metadata['Gene']
    
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
    
    return results

def analyze_feature_importance(results, X, top_n=20):
    """Analyze which motifs are most important for predicting expression"""
    print(f"\nAnalyzing top {top_n} most important motifs...")
    
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
            'Motif': feature_names,
            'Importance': importances,
            'Model': model_name
        }).sort_values('Importance', ascending=False)
        
        importance_data.append(importance_df.head(top_n))
    
    if importance_data:
        all_importance = pd.concat(importance_data, ignore_index=True)
        return all_importance
    else:
        return pd.DataFrame()

def create_visualizations(results, importance_df, merged_df, output_dir):
    """Create visualizations for the analysis"""
    print("Creating visualizations...")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Model performance comparison
    plt.figure(figsize=(10, 6))
    model_names = list(results.keys())
    cv_scores = [results[name]['cv_r2_mean'] for name in model_names]
    cv_errors = [results[name]['cv_r2_std'] for name in model_names]
    
    plt.bar(model_names, cv_scores, yerr=cv_errors, capsize=5)
    plt.ylabel('Cross-Validation R²')
    plt.title('Model Performance Comparison')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/model_performance.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Feature importance heatmap (top motifs across models)
    if not importance_df.empty:
        plt.figure(figsize=(12, 8))
        pivot_importance = importance_df.pivot(index='Motif', columns='Model', values='Importance')
        sns.heatmap(pivot_importance.head(20), annot=True, fmt='.3f', cmap='YlOrRd')
        plt.title('Top 20 Motif Importance Across Models')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/motif_importance_heatmap.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    # 3. Expression distribution by line
    plt.figure(figsize=(10, 6))
    sns.boxplot(data=merged_df, x='Line', y='Expression')
    plt.ylabel('Relative Expression (vs IM767)')
    plt.title('Expression Distribution by Genetic Line')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/expression_by_line.png", dpi=300, bbox_inches='tight')
    plt.close()

def generate_report(results, importance_df, merged_df, output_file):
    """Generate a comprehensive report"""
    print(f"Generating report: {output_file}")
    
    with open(output_file, 'w') as f:
        f.write("# Motif-Expression Regression Analysis Report\n\n")
        
        # Data summary
        f.write("## Data Summary\n")
        f.write(f"- Total gene-line combinations: {len(merged_df)}\n")
        f.write(f"- Unique genes: {merged_df['Gene'].nunique()}\n")
        f.write(f"- Genetic lines: {', '.join(sorted(merged_df['Line'].unique()))}\n")
        f.write(f"- Expression range: {merged_df['Expression'].min():.2f} to {merged_df['Expression'].max():.2f}\n\n")
        
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
        
        # Top motifs
        if not importance_df.empty:
            f.write("## Top 10 Most Important Motifs\n")
            top_motifs = importance_df.groupby('Motif')['Importance'].mean().sort_values(ascending=False).head(10)
            for motif, importance in top_motifs.items():
                clean_motif = motif.replace('_present', '')
                f.write(f"- **{clean_motif}**: {importance:.3f}\n")
            f.write("\n")
        
        # Interpretation
        f.write("## Interpretation\n")
        r2_best = results[best_model]['cv_r2_mean']
        if r2_best > 0.5:
            f.write("🎉 **Strong predictive power**: Motif patterns explain substantial expression variation.\n")
        elif r2_best > 0.2:
            f.write("📈 **Moderate predictive power**: Motif patterns contribute to expression differences.\n")
        else:
            f.write("⚠️ **Limited predictive power**: Consider other regulatory mechanisms or data quality.\n")
        
        f.write(f"\nThe analysis suggests that motif presence patterns can explain "
               f"{r2_best*100:.1f}% of the variance in relative gene expression across lines.\n")

def main():
    parser = argparse.ArgumentParser(
        description='Analyze motif effects on gene expression across genetic lines (ABSOLUTE analysis)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic analysis with binary motif features
  %(prog)s consolidated_motifs.tsv expression_data.tsv --output results/
  
  # Include detailed motif features (counts, positions, etc.)
  %(prog)s consolidated_motifs.tsv expression_data.tsv --detailed --output results/
  
  # Focus on top motifs only
  %(prog)s consolidated_motifs.tsv expression_data.tsv --top-motifs 50 --output results/

Expected file formats:
  consolidated_motifs.tsv: Output from streme_sites_consolidator.py
  expression_data.tsv: 
    Long format: Gene Line Expression (tab-separated)
    Wide format: gene LRTadd IM62 IM155 ... IM767 (tab-separated)

Note: For RELATIVE analysis (comparing to IM767 baseline), use relative_motif_analyzer.py
        """
    )
    
    parser.add_argument('motif_features', help='Consolidated motif file from streme_sites_consolidator.py')
    parser.add_argument('expression_data', help='Expression data file (long or wide format)')
    parser.add_argument('--output', '-o', default='motif_analysis_results', 
                       help='Output directory for results')
    parser.add_argument('--detailed', action='store_true',
                       help='Use detailed motif features (not just presence/absence)')
    parser.add_argument('--top-motifs', type=int, 
                       help='Only use top N most important motifs')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    print("=" * 60)
    print("ABSOLUTE MOTIF-EXPRESSION ANALYSIS")
    print("=" * 60)
    
    try:
        # Load data
        motif_df = load_motif_features(args.motif_features)
        expr_df = load_expression_data(args.expression_data)
        
        # Merge data
        merged_df = merge_motif_expression_data(motif_df, expr_df)
        
        if len(merged_df) == 0:
            print("ERROR: No overlapping data found. Check gene/line names match between files.")
            return 1
        
        # Prepare regression data
        X, y, metadata, motif_cols = prepare_regression_data(merged_df, simple_mode=not args.detailed)
        
        # Filter to top motifs if requested
        if args.top_motifs and args.top_motifs < len(motif_cols):
            # Use variance to select most variable motifs
            motif_variance = X.var().sort_values(ascending=False)
            top_motif_cols = motif_variance.head(args.top_motifs).index
            X = X[top_motif_cols]
            print(f"Filtered to top {args.top_motifs} most variable motifs")
        
        # Run analysis
        results = run_regression_analysis(X, y, metadata)
        
        if not results:
            print("ERROR: No models ran successfully")
            return 1
        
        # Feature importance analysis
        importance_df = analyze_feature_importance(results, X)
        
        # Create visualizations
        create_visualizations(results, importance_df, merged_df, args.output)
        
        # Generate report
        report_file = os.path.join(args.output, 'analysis_report.md')
        generate_report(results, importance_df, merged_df, report_file)
        
        # Save detailed results
        if not importance_df.empty:
            importance_df.to_csv(os.path.join(args.output, 'motif_importance.tsv'), 
                                sep='\t', index=False)
        
        # Save predictions from best model
        best_model_name = max(results.keys(), key=lambda k: results[k]['cv_r2_mean'])
        best_model = results[best_model_name]['model']
        
        predictions_df = metadata.copy()
        predictions_df['Actual_Expression'] = y
        predictions_df['Predicted_Expression'] = best_model.predict(X)
        predictions_df['Residual'] = predictions_df['Actual_Expression'] - predictions_df['Predicted_Expression']
        predictions_df.to_csv(os.path.join(args.output, 'predictions.tsv'), 
                             sep='\t', index=False)
        
        print(f"\n🎉 Analysis complete! Results saved to: {args.output}/")
        print(f"📊 Check analysis_report.md for detailed findings")
        print(f"📈 Best model: {best_model_name} (R² = {results[best_model_name]['cv_r2_mean']:.3f})")
        
        return 0
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())