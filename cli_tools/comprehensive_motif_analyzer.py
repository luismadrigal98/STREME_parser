#!/usr/bin/env python3
"""
Enhanced Motif-Expression Analyzer with Position and Variation Effects

This advanced tool analyzes three layers of regulatory effects:
1. Motif Presence: Basic presence/absence effects
2. Position Bias: How motif position relative to TSS affects expression
3. Sequence Variation: How motif sequence changes across lines affect expression

Key Questions Addressed:
- Do motifs closer to TSS have stronger effects?
- Do sequence variants of the same motif have different regulatory strengths?
- Can we predict expression from the complete regulatory "fingerprint"?
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
import re
import warnings
warnings.filterwarnings('ignore')

def calculate_position_features(motif_sites_df):
    """
    Calculate sophisticated position-based features for each motif in each gene-line
    
    Features:
    - Proximal density (0-500bp from TSS)
    - Core density (500-1000bp)
    - Distal density (1000-2000bp)
    - Position weighted score (closer = higher weight)
    - Spacing regularity (how evenly spaced are motif hits)
    """
    position_features = {}
    
    for (gene, line, motif), group in motif_sites_df.groupby(['gene_id', 'line', 'motif_id']):
        positions = group['position'].values
        
        # Define position bins
        proximal = positions[(positions >= 0) & (positions < 500)]
        core = positions[(positions >= 500) & (positions < 1000)]
        distal = positions[(positions >= 1000) & (positions < 2000)]
        
        # Calculate densities (hits per 100bp)
        proximal_density = len(proximal) / 5.0  # 500bp / 100bp
        core_density = len(core) / 5.0
        distal_density = len(distal) / 10.0    # 1000bp / 100bp
        
        # Position-weighted score (exponential decay from TSS)
        position_weights = np.exp(-positions / 1000.0)  # Decay with 1kb half-life
        weighted_score = np.sum(position_weights)
        
        # Spacing regularity (coefficient of variation of inter-motif distances)
        if len(positions) > 1:
            sorted_pos = np.sort(positions)
            spacings = np.diff(sorted_pos)
            spacing_cv = np.std(spacings) / np.mean(spacings) if np.mean(spacings) > 0 else 0
        else:
            spacing_cv = 0
        
        # Most proximal hit distance
        min_distance = np.min(positions) if len(positions) > 0 else 2000
        
        feature_key = f"{gene}_{line}_{motif}"
        position_features[feature_key] = {
            'proximal_density': proximal_density,
            'core_density': core_density,
            'distal_density': distal_density,
            'position_weighted_score': weighted_score,
            'spacing_regularity': 1.0 / (1.0 + spacing_cv),  # Higher = more regular
            'closest_hit_distance': min_distance,
            'total_hits': len(positions)
        }
    
    return position_features

def calculate_sequence_variation_features(motif_sites_df):
    """
    Calculate sequence variation features for each motif across lines
    
    Features:
    - Consensus deviation: How much each hit deviates from consensus
    - Line-specific consensus: Dominant sequence variant in each line
    - Variation entropy: How diverse are the sequences within each line
    - Cross-line variation: How different are the sequences between lines
    """
    variation_features = {}
    
    # Group by motif to analyze variation patterns
    for motif_id, motif_group in motif_sites_df.groupby('motif_id'):
        consensus = motif_group['consensus'].iloc[0]
        
        # For each gene-line combination
        for (gene, line), gene_line_group in motif_group.groupby(['gene_id', 'line']):
            sequences = gene_line_group['sequence'].values
            
            if len(sequences) == 0:
                continue
            
            # Calculate consensus similarity for each sequence
            similarities = []
            for seq in sequences:
                similarity = calculate_sequence_similarity(seq, consensus)
                similarities.append(similarity)
            
            # Sequence diversity within this gene-line
            unique_seqs = len(set(sequences))
            diversity = unique_seqs / len(sequences) if len(sequences) > 0 else 0
            
            # GC content variation
            gc_contents = [calculate_gc_content(seq) for seq in sequences]
            gc_mean = np.mean(gc_contents)
            gc_std = np.std(gc_contents)
            
            # Sequence length variation
            lengths = [len(seq) for seq in sequences]
            length_mean = np.mean(lengths)
            length_std = np.std(lengths)
            
            feature_key = f"{gene}_{line}_{motif_id}"
            variation_features[feature_key] = {
                'consensus_similarity_mean': np.mean(similarities),
                'consensus_similarity_std': np.std(similarities),
                'sequence_diversity': diversity,
                'gc_content_mean': gc_mean,
                'gc_content_std': gc_std,
                'length_mean': length_mean,
                'length_std': length_std,
                'unique_sequences': unique_seqs
            }
    
    return variation_features

def calculate_cross_line_variation_features(motif_sites_df):
    """
    Calculate how motif sequences vary across genetic lines for each gene-motif combination
    """
    cross_line_features = {}
    
    # Group by gene and motif to compare across lines
    for (gene, motif_id), gene_motif_group in motif_sites_df.groupby(['gene_id', 'motif_id']):
        
        # Collect sequences from each line
        line_sequences = {}
        for line, line_group in gene_motif_group.groupby('line'):
            line_sequences[line] = line_group['sequence'].tolist()
        
        # Calculate cross-line variation
        all_sequences = []
        for line_seqs in line_sequences.values():
            all_sequences.extend(line_seqs)
        
        if len(all_sequences) == 0:
            continue
        
        # Overall diversity across all lines
        total_unique = len(set(all_sequences))
        total_sequences = len(all_sequences)
        cross_line_diversity = total_unique / total_sequences if total_sequences > 0 else 0
        
        # Line-specific dominant sequences
        line_dominant_seqs = {}
        for line, seqs in line_sequences.items():
            if seqs:
                # Find most common sequence in this line
                seq_counts = pd.Series(seqs).value_counts()
                line_dominant_seqs[line] = seq_counts.index[0]
        
        # Calculate how different the dominant sequences are between lines
        dominant_seqs = list(line_dominant_seqs.values())
        if len(dominant_seqs) > 1:
            dominant_diversity = len(set(dominant_seqs)) / len(dominant_seqs)
        else:
            dominant_diversity = 0
        
        # For each line, calculate features
        for line in line_sequences.keys():
            feature_key = f"{gene}_{line}_{motif_id}"
            cross_line_features[feature_key] = {
                'cross_line_diversity': cross_line_diversity,
                'dominant_sequence_diversity': dominant_diversity,
                'line_has_unique_variant': int(line in line_dominant_seqs and 
                                             line_dominant_seqs[line] not in 
                                             [seq for other_line, seq in line_dominant_seqs.items() 
                                              if other_line != line])
            }
    
    return cross_line_features

def calculate_sequence_similarity(seq1, seq2):
    """Calculate similarity between two sequences with IUPAC awareness"""
    if len(seq1) != len(seq2):
        return 0.0
    
    matches = 0
    for i, (base1, base2) in enumerate(zip(seq1.upper(), seq2.upper())):
        if base1 == base2 or base1 == 'N' or base2 == 'N':
            matches += 1
        elif base1 in 'RYSWKMBDHV' or base2 in 'RYSWKMBDHV':
            # IUPAC code matching - simplified
            matches += 0.5
    
    return matches / len(seq1)

def calculate_gc_content(sequence):
    """Calculate GC content of a sequence"""
    if not sequence:
        return 0.0
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    return gc_count / len(sequence)

def load_consolidated_motif_data(motif_file):
    """Load consolidated STREME sites data"""
    print(f"Loading consolidated motif data from: {motif_file}")
    
    df = pd.read_csv(motif_file, sep='\t')
    
    # Ensure required columns exist
    required_cols = ['gene_id', 'line', 'motif_id', 'consensus', 'sequence', 'position']
    missing_cols = [col for col in required_cols if col not in df.columns]
    if missing_cols:
        print(f"ERROR: Missing required columns: {missing_cols}")
        return None
    
    print(f"Loaded {len(df)} motif sites")
    print(f"Genes: {df['gene_id'].nunique()}")
    print(f"Lines: {df['line'].nunique()}")
    print(f"Motifs: {df['motif_id'].nunique()}")
    
    return df

def create_comprehensive_features(motif_sites_df, feature_types=['presence', 'position', 'variation', 'cross_line']):
    """
    Create comprehensive feature matrix including presence, position, and variation features
    """
    print("Creating comprehensive feature matrix...")
    
    # Initialize feature collections
    all_features = defaultdict(dict)
    
    # 1. Basic presence features
    if 'presence' in feature_types:
        print("  - Adding presence features...")
        for (gene, line, motif), group in motif_sites_df.groupby(['gene_id', 'line', 'motif_id']):
            key = f"{gene}_{line}"
            all_features[key][f"{motif}_present"] = 1
            all_features[key][f"{motif}_count"] = len(group)
    
    # 2. Position-based features
    if 'position' in feature_types:
        print("  - Adding position features...")
        position_features = calculate_position_features(motif_sites_df)
        for feature_key, features in position_features.items():
            gene, line, motif = feature_key.split('_', 2)
            key = f"{gene}_{line}"
            for feat_name, feat_value in features.items():
                all_features[key][f"{motif}_{feat_name}"] = feat_value
    
    # 3. Sequence variation features
    if 'variation' in feature_types:
        print("  - Adding sequence variation features...")
        variation_features = calculate_sequence_variation_features(motif_sites_df)
        for feature_key, features in variation_features.items():
            gene, line, motif = feature_key.split('_', 2)
            key = f"{gene}_{line}"
            for feat_name, feat_value in features.items():
                all_features[key][f"{motif}_{feat_name}"] = feat_value
    
    # 4. Cross-line variation features
    if 'cross_line' in feature_types:
        print("  - Adding cross-line variation features...")
        cross_line_features = calculate_cross_line_variation_features(motif_sites_df)
        for feature_key, features in cross_line_features.items():
            gene, line, motif = feature_key.split('_', 2)
            key = f"{gene}_{line}"
            for feat_name, feat_value in features.items():
                all_features[key][f"{motif}_{feat_name}"] = feat_value
    
    # Convert to DataFrame
    feature_df = pd.DataFrame.from_dict(all_features, orient='index')
    feature_df = feature_df.fillna(0)  # Fill missing values with 0
    
    # Extract gene and line from index
    feature_df['Gene'] = feature_df.index.str.split('_').str[0]
    feature_df['Line'] = feature_df.index.str.split('_').str[1]
    
    print(f"Created feature matrix: {feature_df.shape[0]} samples × {feature_df.shape[1]-2} features")
    
    return feature_df

def analyze_regulatory_layers(X, y, metadata, feature_types):
    """
    Analyze the contribution of different regulatory layers
    """
    print("\nAnalyzing regulatory layer contributions...")
    
    results = {}
    
    # Analyze each layer separately
    for layer in feature_types:
        layer_cols = [col for col in X.columns if any(suffix in col for suffix in get_layer_suffixes(layer))]
        
        if not layer_cols:
            continue
        
        X_layer = X[layer_cols]
        
        # Use a simple model for layer analysis
        model = RandomForestRegressor(n_estimators=100, random_state=42)
        cv = GroupKFold(n_splits=5)
        groups = metadata['Gene']
        
        cv_scores = cross_val_score(model, X_layer, y, cv=cv, groups=groups, scoring='r2')
        
        results[layer] = {
            'r2_mean': cv_scores.mean(),
            'r2_std': cv_scores.std(),
            'n_features': len(layer_cols)
        }
        
        print(f"{layer:15s}: R² = {cv_scores.mean():.3f} ± {cv_scores.std():.3f} ({len(layer_cols)} features)")
    
    return results

def get_layer_suffixes(layer):
    """Get feature suffixes for each regulatory layer"""
    suffixes = {
        'presence': ['_present', '_count'],
        'position': ['_proximal_density', '_core_density', '_distal_density', 
                    '_position_weighted_score', '_spacing_regularity', '_closest_hit_distance'],
        'variation': ['_consensus_similarity_mean', '_consensus_similarity_std', 
                     '_sequence_diversity', '_gc_content_mean', '_gc_content_std'],
        'cross_line': ['_cross_line_diversity', '_dominant_sequence_diversity', 
                      '_line_has_unique_variant']
    }
    return suffixes.get(layer, [])

def create_interaction_features(X, top_motifs=10):
    """
    Create interaction features between motifs, positions, and variations
    """
    print(f"Creating interaction features for top {top_motifs} motifs...")
    
    # Get top motifs by variance
    presence_cols = [col for col in X.columns if col.endswith('_present')]
    if len(presence_cols) == 0:
        return X
    
    motif_variance = X[presence_cols].var().sort_values(ascending=False)
    top_motif_ids = [col.replace('_present', '') for col in motif_variance.head(top_motifs).index]
    
    interaction_features = {}
    
    for i, motif1 in enumerate(top_motif_ids):
        for j, motif2 in enumerate(top_motif_ids[i+1:], i+1):
            # Presence * Presence interaction
            if f"{motif1}_present" in X.columns and f"{motif2}_present" in X.columns:
                interaction_features[f"{motif1}_X_{motif2}_presence"] = (
                    X[f"{motif1}_present"] * X[f"{motif2}_present"]
                )
            
            # Position * Presence interaction
            if f"{motif1}_position_weighted_score" in X.columns and f"{motif2}_present" in X.columns:
                interaction_features[f"{motif1}_position_X_{motif2}_presence"] = (
                    X[f"{motif1}_position_weighted_score"] * X[f"{motif2}_present"]
                )
    
    # Add to original features
    interaction_df = pd.DataFrame(interaction_features, index=X.index)
    X_enhanced = pd.concat([X, interaction_df], axis=1)
    
    print(f"Added {len(interaction_features)} interaction features")
    
    return X_enhanced

def main():
    parser = argparse.ArgumentParser(
        description='Comprehensive motif-expression analysis with position and variation effects',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Full analysis with all regulatory layers
  %(prog)s consolidated_sites.tsv expression_data.tsv --output results/
  
  # Focus on position effects only
  %(prog)s consolidated_sites.tsv expression_data.tsv --layers presence position --output results/
  
  # Include interaction effects
  %(prog)s consolidated_sites.tsv expression_data.tsv --interactions --output results/
        """
    )
    
    parser.add_argument('motif_sites_file', help='Consolidated STREME sites TSV file')
    parser.add_argument('expression_data', help='Expression data file (Gene, Line, Expression)')
    parser.add_argument('--output', '-o', default='comprehensive_analysis', 
                       help='Output directory for results')
    parser.add_argument('--layers', nargs='+', 
                       choices=['presence', 'position', 'variation', 'cross_line'],
                       default=['presence', 'position', 'variation', 'cross_line'],
                       help='Regulatory layers to analyze')
    parser.add_argument('--interactions', action='store_true',
                       help='Include motif-motif and motif-position interactions')
    parser.add_argument('--top-motifs', type=int, default=50,
                       help='Number of top motifs to focus analysis on')
    
    args = parser.parse_args()
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    print("=" * 70)
    print("COMPREHENSIVE MOTIF-EXPRESSION ANALYSIS")
    print("Position Bias + Sequence Variation + Cross-Line Effects")
    print("=" * 70)
    
    try:
        # Load motif sites data
        motif_sites_df = load_consolidated_motif_data(args.motif_sites_file)
        if motif_sites_df is None:
            return 1
        
        # Load expression data
        expr_df = pd.read_csv(args.expression_data, sep='\t')
        expected_cols = ['Gene', 'Line', 'Expression']
        if list(expr_df.columns) != expected_cols:
            expr_df.columns = expected_cols[:len(expr_df.columns)]
        
        # Create comprehensive features
        feature_df = create_comprehensive_features(motif_sites_df, args.layers)
        
        # Merge with expression data
        merged_df = pd.merge(feature_df, expr_df, on=['Gene', 'Line'], how='inner')
        
        if len(merged_df) == 0:
            print("ERROR: No overlapping data found.")
            return 1
        
        print(f"Merged data: {len(merged_df)} gene-line combinations")
        
        # Prepare for analysis
        feature_cols = [col for col in merged_df.columns if col not in ['Gene', 'Line', 'Expression']]
        X = merged_df[feature_cols]
        y = merged_df['Expression']
        metadata = merged_df[['Gene', 'Line']]
        
        # Filter to top variable motifs
        if args.top_motifs:
            # Get most variable motifs
            presence_cols = [col for col in X.columns if col.endswith('_present')]
            if presence_cols:
                motif_variance = X[presence_cols].var().sort_values(ascending=False)
                top_motifs = [col.replace('_present', '') for col in motif_variance.head(args.top_motifs).index]
                
                # Keep only features for top motifs
                top_feature_cols = [col for col in X.columns if any(motif in col for motif in top_motifs)]
                X = X[top_feature_cols]
                print(f"Filtered to top {args.top_motifs} motifs: {X.shape[1]} features")
        
        # Add interaction features if requested
        if args.interactions:
            X = create_interaction_features(X, top_motifs=min(10, args.top_motifs or 10))
        
        # Analyze regulatory layer contributions
        layer_results = analyze_regulatory_layers(X, y, metadata, args.layers)
        
        # Run comprehensive model
        print("\nRunning comprehensive model...")
        model = RandomForestRegressor(n_estimators=200, random_state=42)
        cv = GroupKFold(n_splits=5)
        groups = metadata['Gene']
        
        cv_scores = cross_val_score(model, X, y, cv=cv, groups=groups, scoring='r2')
        print(f"Comprehensive model: R² = {cv_scores.mean():.3f} ± {cv_scores.std():.3f}")
        
        # Feature importance analysis
        model.fit(X, y)
        feature_importance = pd.DataFrame({
            'Feature': X.columns,
            'Importance': model.feature_importances_
        }).sort_values('Importance', ascending=False)
        
        # Save results
        feature_importance.to_csv(f"{args.output}/feature_importance.tsv", sep='\t', index=False)
        
        # Create summary report
        with open(f"{args.output}/analysis_summary.md", 'w') as f:
            f.write("# Comprehensive Motif-Expression Analysis Results\n\n")
            
            f.write("## Regulatory Layer Contributions\n")
            for layer, results in layer_results.items():
                f.write(f"- **{layer.title()}**: R² = {results['r2_mean']:.3f} ± {results['r2_std']:.3f} "
                       f"({results['n_features']} features)\n")
            
            f.write(f"\n## Comprehensive Model Performance\n")
            f.write(f"**Combined Model**: R² = {cv_scores.mean():.3f} ± {cv_scores.std():.3f}\n\n")
            
            f.write("## Top 20 Most Important Features\n")
            for _, row in feature_importance.head(20).iterrows():
                feature_type = "Unknown"
                if "_present" in row['Feature']:
                    feature_type = "Presence"
                elif any(pos in row['Feature'] for pos in ['proximal', 'core', 'distal', 'position']):
                    feature_type = "Position"
                elif any(var in row['Feature'] for var in ['similarity', 'diversity', 'gc_content']):
                    feature_type = "Variation"
                elif "cross_line" in row['Feature']:
                    feature_type = "Cross-line"
                elif "_X_" in row['Feature']:
                    feature_type = "Interaction"
                
                f.write(f"- **{row['Feature']}** ({feature_type}): {row['Importance']:.4f}\n")
        
        print(f"\n🎉 Analysis complete! Results saved to: {args.output}/")
        print(f"📊 Check analysis_summary.md for detailed findings")
        print(f"📈 Comprehensive model R² = {cv_scores.mean():.3f}")
        
        return 0
        
    except Exception as e:
        print(f"ERROR: {e}")
        import traceback
        traceback.print_exc()
        return 1

if __name__ == "__main__":
    sys.exit(main())