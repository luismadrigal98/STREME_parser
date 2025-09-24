#!/usr/bin/env python3
"""
Gene-Specific Motif-Expression Analysis Tool

This tool performs detailed motif-expression analysis for individual genes or gene sets,
using all available motifs in their promoter regions. This approach avoids the 
high-dimensional problem of genome-wide analysis while providing detailed insights
into gene-specific regulatory mechanisms.

Features:
- Focus on individual genes or small gene sets
- Use all motifs present in target gene promoters
- Support both absolute and relative analysis
- Detailed motif-by-motif contribution analysis
- Reduced multiple testing burden

Usage:
    python gene_specific_analyzer.py <consolidated_motifs.tsv> <expression_data.tsv> --genes <gene_list> [options]
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
from sklearn.model_selection import cross_val_score, LeaveOneOut
from sklearn.metrics import r2_score, mean_squared_error
from sklearn.preprocessing import StandardScaler
from collections import defaultdict
import warnings
warnings.filterwarnings('ignore')

def load_motif_data(motif_file):
    """Load consolidated motif data"""
    print(f"Loading consolidated motif data from: {motif_file}")
    
    df = pd.read_csv(motif_file, sep='\t')
    
    print(f"Loaded {len(df)} motif sites")
    print(f"Genes: {df['gene_id'].nunique()}")
    print(f"Lines: {df['line'].nunique()}")
    print(f"Motifs: {df['consolidated_motif_id'].nunique()}")
    
    return df

def load_expression_data(expression_file):
    """Load expression data in wide or long format"""
    print(f"Loading expression data from: {expression_file}")
    
    expr_df = pd.read_csv(expression_file, sep='\t')
    
    # Detect and convert format if needed (same logic as before)
    if len(expr_df.columns) > 3 and any(col.startswith('IM') for col in expr_df.columns):
        print("Detected wide format - converting to long format")
        
        gene_col = expr_df.columns[0]
        line_cols = [col for col in expr_df.columns if col not in [gene_col, 'LRTadd'] and col.startswith('IM')]
        
        long_df = expr_df.melt(
            id_vars=[gene_col], 
            value_vars=line_cols,
            var_name='Line', 
            value_name='Expression'
        )
        long_df = long_df.rename(columns={gene_col: 'Gene'})
        
        print(f"Found expression data for {long_df['Gene'].nunique()} genes across {long_df['Line'].nunique()} lines")
        return long_df
    
    else:
        # Long format
        expected_cols = ['Gene', 'Line', 'Expression']
        if list(expr_df.columns) != expected_cols:
            print(f"Assuming first 3 columns are Gene, Line, Expression")
            expr_df.columns = expected_cols[:len(expr_df.columns)]
        
        print(f"Found expression data for {expr_df['Gene'].nunique()} genes across {expr_df['Line'].nunique()} lines")
        return expr_df

def extract_gene_motif_features(motif_df, target_genes, analysis_type='relative', reference_line='IM767'):
    """
    Extract comprehensive motif features for target genes
    """
    print(f"Extracting motif features for {len(target_genes)} target genes...")
    
    # Filter to target genes
    gene_motifs = motif_df[motif_df['gene_id'].isin(target_genes)].copy()
    
    if len(gene_motifs) == 0:
        print("WARNING: No motif data found for target genes")
        return pd.DataFrame()
    
    print(f"Found motif data for {gene_motifs['gene_id'].nunique()} of {len(target_genes)} target genes")
    print(f"Total motif sites: {len(gene_motifs)}")
    print(f"Unique motifs: {gene_motifs['consolidated_motif_id'].nunique()}")
    
    # Calculate comprehensive features per gene-line-motif
    features_list = []
    
    for gene in gene_motifs['gene_id'].unique():
        gene_data = gene_motifs[gene_motifs['gene_id'] == gene]
        
        for line in gene_data['line'].unique():
            line_data = gene_data[gene_data['line'] == line]
            
            # Get all motifs for this gene-line combo
            motif_features = {'gene': gene, 'line': line}
            
            for motif in gene_data['consolidated_motif_id'].unique():
                motif_sites = line_data[line_data['consolidated_motif_id'] == motif]
                
                prefix = f"motif_{motif}"
                
                if len(motif_sites) > 0:
                    # Presence/count features
                    motif_features[f"{prefix}_present"] = 1
                    motif_features[f"{prefix}_count"] = len(motif_sites)
                    
                    # Position features
                    positions = motif_sites['start_pos'].values
                    motif_features[f"{prefix}_position_mean"] = np.mean(positions)
                    motif_features[f"{prefix}_position_std"] = np.std(positions) if len(positions) > 1 else 0
                    motif_features[f"{prefix}_position_min"] = np.min(positions)
                    motif_features[f"{prefix}_position_max"] = np.max(positions)
                    
                    # Score features
                    if 'score' in motif_sites.columns:
                        scores = motif_sites['score'].values
                        motif_features[f"{prefix}_score_mean"] = np.mean(scores)
                        motif_features[f"{prefix}_score_max"] = np.max(scores)
                    
                    # Sequence features
                    if 'sequence' in motif_sites.columns:
                        sequences = motif_sites['sequence'].values
                        gc_contents = [(seq.count('G') + seq.count('C')) / len(seq) for seq in sequences if isinstance(seq, str)]
                        if gc_contents:
                            motif_features[f"{prefix}_gc_content"] = np.mean(gc_contents)
                
                else:
                    # Absent motif
                    motif_features[f"{prefix}_present"] = 0
                    motif_features[f"{prefix}_count"] = 0
                    motif_features[f"{prefix}_position_mean"] = 0
                    motif_features[f"{prefix}_position_std"] = 0
                    motif_features[f"{prefix}_position_min"] = 0  
                    motif_features[f"{prefix}_position_max"] = 0
                    motif_features[f"{prefix}_score_mean"] = 0
                    motif_features[f"{prefix}_score_max"] = 0
                    motif_features[f"{prefix}_gc_content"] = 0
            
            features_list.append(motif_features)
    
    features_df = pd.DataFrame(features_list)
    
    # Convert to relative features if requested
    if analysis_type == 'relative':
        print(f"Converting to relative features (reference: {reference_line})")
        relative_features = []
        
        for gene in features_df['gene'].unique():
            gene_data = features_df[features_df['gene'] == gene]
            ref_data = gene_data[gene_data['line'] == reference_line]
            
            if len(ref_data) == 0:
                print(f"WARNING: No reference line data for gene {gene}")
                continue
            
            ref_row = ref_data.iloc[0]
            
            for _, comp_row in gene_data[gene_data['line'] != reference_line].iterrows():
                relative_row = {'gene': gene, 'comparison_line': comp_row['line']}
                
                # Calculate differences for all motif features
                for col in gene_data.columns:
                    if col.startswith('motif_'):
                        relative_row[f"{col}_diff"] = comp_row[col] - ref_row[col]
                
                relative_features.append(relative_row)
        
        features_df = pd.DataFrame(relative_features)
    
    print(f"Generated feature matrix: {features_df.shape}")
    return features_df

def run_gene_specific_analysis(features_df, expr_df, analysis_type='relative', reference_line='IM767'):
    """
    Run regression analysis on gene-specific features
    """
    print(f"\nRunning {analysis_type} gene-specific analysis...")
    
    # Merge with expression data
    if analysis_type == 'relative':
        # For relative analysis, create expression differences
        expr_diffs = []
        for gene in features_df['gene'].unique():
            gene_expr = expr_df[expr_df['Gene'] == gene]
            ref_expr = gene_expr[gene_expr['Line'] == reference_line]
            
            if len(ref_expr) == 0:
                continue
            
            ref_value = ref_expr['Expression'].iloc[0]
            
            for _, row in gene_expr[gene_expr['Line'] != reference_line].iterrows():
                expr_diffs.append({
                    'gene': gene,
                    'comparison_line': row['Line'],
                    'expression_difference': row['Expression'] - ref_value
                })
        
        expr_df_processed = pd.DataFrame(expr_diffs)
        merged_df = pd.merge(features_df, expr_df_processed, on=['gene', 'comparison_line'], how='inner')
        target_col = 'expression_difference'
        
    else:
        # Absolute analysis
        expr_df_processed = expr_df.rename(columns={'Gene': 'gene', 'Line': 'line'})
        merged_df = pd.merge(features_df, expr_df_processed, on=['gene', 'line'], how='inner')
        target_col = 'Expression'
    
    if len(merged_df) == 0:
        print("ERROR: No overlapping data between motif features and expression")
        return {}
    
    print(f"Analyzing {len(merged_df)} gene-line combinations")
    
    # Prepare feature matrix
    feature_cols = [col for col in merged_df.columns if col.startswith('motif_')]
    X = merged_df[feature_cols].copy()
    y = merged_df[target_col].copy()
    
    # Remove constant features
    constant_features = X.columns[X.var() == 0]
    if len(constant_features) > 0:
        print(f"Removing {len(constant_features)} constant features")
        X = X.drop(columns=constant_features)
    
    print(f"Final feature matrix: {X.shape[0]} samples × {X.shape[1]} features")
    print(f"Expression range: {y.min():.3f} to {y.max():.3f}")
    
    if X.shape[1] == 0:
        print("ERROR: No variable features found")
        return {}
    
    # Run models with appropriate cross-validation
    if len(X) <= 10:
        # Use Leave-One-Out for small datasets
        cv = LeaveOneOut()
        cv_name = "Leave-One-Out"
    else:
        # Use 5-fold for larger datasets
        from sklearn.model_selection import KFold
        cv = KFold(n_splits=min(5, len(X)), shuffle=True, random_state=42)
        cv_name = f"{min(5, len(X))}-Fold"
    
    models = {
        'Linear Regression': LinearRegression(),
        'Ridge Regression': RidgeCV(alphas=[0.1, 1.0, 10.0, 100.0]),
        'Lasso Regression': LassoCV(cv=3, max_iter=10000),
        'Random Forest': RandomForestRegressor(n_estimators=100, random_state=42),
    }
    
    results = {}
    
    print(f"\nUsing {cv_name} cross-validation:")
    
    for name, model in models.items():
        try:
            # Cross-validation scores
            cv_scores = cross_val_score(model, X, y, cv=cv, scoring='r2')
            
            # Fit full model
            model.fit(X, y)
            y_pred = model.predict(X)
            
            results[name] = {
                'cv_r2_mean': cv_scores.mean(),
                'cv_r2_std': cv_scores.std(),
                'full_r2': r2_score(y, y_pred),
                'rmse': np.sqrt(mean_squared_error(y, y_pred)),
                'model': model,
                'feature_names': X.columns.tolist()
            }
            
            print(f"{name:20s} - CV R²: {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
            
        except Exception as e:
            print(f"Error with {name}: {e}")
            continue
    
    return results, merged_df, X, y

def create_gene_specific_visualizations(results, merged_df, X, y, output_dir, analysis_type):
    """Create visualizations for gene-specific analysis"""
    print("Creating gene-specific visualizations...")
    
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Model performance
    plt.figure(figsize=(10, 6))
    model_names = list(results.keys())
    cv_scores = [results[name]['cv_r2_mean'] for name in model_names]
    cv_errors = [results[name]['cv_r2_std'] for name in model_names]
    
    plt.bar(model_names, cv_scores, yerr=cv_errors, capsize=5)
    plt.ylabel('Cross-Validation R²')
    plt.title(f'Gene-Specific {analysis_type.title()} Analysis - Model Performance')
    plt.xticks(rotation=45)
    plt.tight_layout()
    plt.savefig(f"{output_dir}/model_performance.png", dpi=300, bbox_inches='tight')
    plt.close()
    
    # 2. Feature importance (if available)
    best_model_name = max(results.keys(), key=lambda k: results[k]['cv_r2_mean'])
    best_model = results[best_model_name]['model']
    
    if hasattr(best_model, 'feature_importances_'):
        importance = best_model.feature_importances_
    elif hasattr(best_model, 'coef_'):
        importance = np.abs(best_model.coef_)
    else:
        importance = None
    
    if importance is not None:
        feature_names = X.columns
        importance_df = pd.DataFrame({
            'Feature': feature_names,
            'Importance': importance
        }).sort_values('Importance', ascending=True)
        
        # Plot top features
        top_features = importance_df.tail(min(20, len(importance_df)))
        
        plt.figure(figsize=(10, max(6, len(top_features) * 0.3)))
        plt.barh(top_features['Feature'], top_features['Importance'])
        plt.xlabel('Feature Importance')
        plt.title(f'Top Features - {best_model_name}')
        plt.tight_layout()
        plt.savefig(f"{output_dir}/feature_importance.png", dpi=300, bbox_inches='tight')
        plt.close()
    
    # 3. Predictions vs Actual
    y_pred = best_model.predict(X)
    plt.figure(figsize=(8, 8))
    plt.scatter(y, y_pred, alpha=0.6)
    plt.plot([y.min(), y.max()], [y.min(), y.max()], 'r--', lw=2)
    plt.xlabel('Actual Expression')
    plt.ylabel('Predicted Expression')
    plt.title(f'Predictions vs Actual - {best_model_name}\nR² = {results[best_model_name]["full_r2"]:.3f}')
    plt.tight_layout()
    plt.savefig(f"{output_dir}/predictions_vs_actual.png", dpi=300, bbox_inches='tight')
    plt.close()

def generate_gene_specific_report(results, merged_df, target_genes, output_file, analysis_type, reference_line):
    """Generate comprehensive report for gene-specific analysis"""
    print(f"Generating gene-specific report: {output_file}")
    
    with open(output_file, 'w') as f:
        f.write(f"# Gene-Specific {analysis_type.title()} Motif-Expression Analysis Report\n\n")
        
        # Target genes
        f.write("## Target Genes\n")
        f.write(f"Analyzed genes: {', '.join(target_genes)}\n\n")
        
        if analysis_type == 'relative':
            f.write(f"Reference line: {reference_line}\n\n")
        
        # Data summary
        f.write("## Data Summary\n")
        f.write(f"- Total observations: {len(merged_df)}\n")
        f.write(f"- Unique genes with data: {merged_df['gene'].nunique()}\n")
        
        if analysis_type == 'relative':
            f.write(f"- Comparison lines: {merged_df['comparison_line'].nunique()}\n")
        else:
            f.write(f"- Genetic lines: {merged_df['line'].nunique()}\n")
        
        # Model performance
        f.write("\n## Model Performance\n")
        f.write("| Model | CV R² | CV Std | Full R² | RMSE |\n")
        f.write("|-------|-------|--------|---------|------|\n")
        for name, result in results.items():
            f.write(f"| {name} | {result['cv_r2_mean']:.3f} | {result['cv_r2_std']:.3f} | "
                   f"{result['full_r2']:.3f} | {result['rmse']:.3f} |\n")
        
        # Best model
        best_model_name = max(results.keys(), key=lambda k: results[k]['cv_r2_mean'])
        f.write(f"\n**Best performing model:** {best_model_name} (CV R² = {results[best_model_name]['cv_r2_mean']:.3f})\n\n")
        
        # Interpretation
        f.write("## Interpretation\n")
        r2_best = results[best_model_name]['cv_r2_mean']
        if r2_best > 0.7:
            f.write("🎉 **Excellent predictive power**: Motif patterns strongly predict expression in target genes.\n")
        elif r2_best > 0.4:
            f.write("📈 **Good predictive power**: Motif patterns substantially contribute to expression variation.\n")
        elif r2_best > 0.1:
            f.write("📊 **Moderate predictive power**: Some motif-expression relationships detected.\n")
        else:
            f.write("⚠️ **Limited predictive power**: Consider data quality or alternative regulatory mechanisms.\n")
        
        f.write(f"\nGene-specific analysis suggests that motif patterns can explain "
               f"{r2_best*100:.1f}% of expression variation in the target genes.\n")

def main():
    parser = argparse.ArgumentParser(
        description='Gene-specific motif-expression analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Analyze specific genes with relative approach
  %(prog)s consolidated_motifs.tsv expression.tsv --genes AT1G01010,AT1G01020 --type relative
  
  # Analyze genes from file with absolute approach  
  %(prog)s consolidated_motifs.tsv expression.tsv --gene-file target_genes.txt --type absolute
  
  # Focus on highly expressed genes
  %(prog)s consolidated_motifs.tsv expression.tsv --genes AT1G01010 --type relative --reference-line IM500
        """
    )
    
    parser.add_argument('consolidated_file', help='Consolidated motif TSV file')
    parser.add_argument('expression_file', help='Expression data file')
    parser.add_argument('--genes', help='Comma-separated list of target genes')
    parser.add_argument('--gene-file', help='File containing one gene ID per line')
    parser.add_argument('--type', choices=['absolute', 'relative'], default='relative',
                       help='Analysis type (default: relative)')
    parser.add_argument('--reference-line', '-r', default='IM767',
                       help='Reference line for relative analysis (default: IM767)')
    parser.add_argument('--output', '-o', default='gene_specific_results/',
                       help='Output directory (default: gene_specific_results/)')
    
    args = parser.parse_args()
    
    # Get target genes
    if args.genes:
        target_genes = [gene.strip() for gene in args.genes.split(',')]
    elif args.gene_file:
        with open(args.gene_file) as f:
            target_genes = [line.strip() for line in f if line.strip()]
    else:
        print("ERROR: Must specify either --genes or --gene-file")
        return 1
    
    print(f"Target genes ({len(target_genes)}): {', '.join(target_genes[:5])}{'...' if len(target_genes) > 5 else ''}")
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    try:
        # Load data
        motif_df = load_motif_data(args.consolidated_file)
        expr_df = load_expression_data(args.expression_file)
        
        # Extract gene-specific features
        features_df = extract_gene_motif_features(
            motif_df, target_genes, args.type, args.reference_line
        )
        
        if len(features_df) == 0:
            print("ERROR: No features generated for target genes")
            return 1
        
        # Run analysis
        results, merged_df, X, y = run_gene_specific_analysis(
            features_df, expr_df, args.type, args.reference_line
        )
        
        if not results:
            print("ERROR: No models ran successfully")
            return 1
        
        # Create visualizations
        create_gene_specific_visualizations(
            results, merged_df, X, y, args.output, args.type
        )
        
        # Generate report
        report_file = os.path.join(args.output, 'gene_specific_analysis_report.md')
        generate_gene_specific_report(
            results, merged_df, target_genes, report_file, args.type, args.reference_line
        )
        
        # Save detailed results
        features_df.to_csv(os.path.join(args.output, 'gene_specific_features.tsv'), 
                          sep='\t', index=False)
        merged_df.to_csv(os.path.join(args.output, 'analysis_data.tsv'), 
                        sep='\t', index=False)
        
        print(f"\n🎉 Gene-specific analysis complete!")
        print(f"📊 Results saved to: {args.output}/")
        best_model_name = max(results.keys(), key=lambda k: results[k]['cv_r2_mean'])
        print(f"📈 Best model: {best_model_name} (R² = {results[best_model_name]['cv_r2_mean']:.3f})")
        
        return 0
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())