#!/usr/bin/env python3
"""
Convert consolidated motif data into regression-ready features for gene expression prediction.
Creates features that capture motif presence, sequence variation, and positional information.
"""

import csv
import argparse
import sys
import os
from collections import defaultdict, Counter
import numpy as np
from difflib import SequenceMatcher

def calculate_sequence_similarity(seq1, seq2):
    """Calculate similarity between two sequences"""
    if not seq1 or not seq2:
        return 0.0
    return SequenceMatcher(None, seq1.upper(), seq2.upper()).ratio()

def calculate_position_features(positions, promoter_length=1000):
    """
    Calculate positional features for motif occurrences
    
    Args:
        positions: List of motif positions
        promoter_length: Length of promoter region (default 1000bp upstream)
    
    Returns:
        Dictionary of positional features
    """
    if not positions:
        return {
            'motif_count': 0,
            'position_mean': 0,
            'position_std': 0,
            'position_min': 0,
            'position_max': 0,
            'position_range': 0,
            'density_proximal': 0,  # Close to TSS (700-1000bp)
            'density_core': 0,      # Middle region (300-700bp)
            'density_distal': 0,    # Far from TSS (0-300bp)
            'spacing_regularity': 0 # How regularly spaced are motifs
        }
    
    positions = [abs(int(p)) for p in positions]  # Convert to absolute positions
    
    # Basic statistics
    pos_mean = np.mean(positions)
    pos_std = np.std(positions) if len(positions) > 1 else 0
    pos_min = min(positions)
    pos_max = max(positions)
    pos_range = pos_max - pos_min
    
    # Regional density (motifs per 100bp)
    # For 1000bp upstream: proximal = close to TSS (700-1000), distal = far from TSS (0-300)
    proximal_count = sum(1 for p in positions if 700 <= p <= 1000)  # Close to gene/TSS
    core_count = sum(1 for p in positions if 300 <= p < 700)        # Middle region
    distal_count = sum(1 for p in positions if 0 <= p < 300)        # Far from gene
    
    density_proximal = proximal_count / 3.0  # Per 100bp (300bp region)
    density_core = core_count / 4.0          # Per 100bp (400bp region)  
    density_distal = distal_count / 3.0      # Per 100bp (300bp region)
    
    # Spacing regularity (lower values = more regular spacing)
    spacing_regularity = 0
    if len(positions) > 2:
        sorted_pos = sorted(positions)
        spacings = [sorted_pos[i+1] - sorted_pos[i] for i in range(len(sorted_pos)-1)]
        spacing_regularity = np.std(spacings) / (np.mean(spacings) + 1)  # Coefficient of variation
    
    return {
        'motif_count': len(positions),
        'position_mean': pos_mean,
        'position_std': pos_std,
        'position_min': pos_min,
        'position_max': pos_max,
        'position_range': pos_range,
        'density_proximal': density_proximal,
        'density_core': density_core,
        'density_distal': density_distal,
        'spacing_regularity': spacing_regularity
    }

def calculate_sequence_variation_features(sequences, consensus):
    """
    Calculate features that capture sequence variation within a motif family
    
    Args:
        sequences: List of actual motif sequences
        consensus: Consensus sequence for this motif
    
    Returns:
        Dictionary of variation features
    """
    if not sequences:
        return {
            'seq_diversity': 0,
            'consensus_similarity_mean': 0,
            'consensus_similarity_std': 0,
            'seq_length_mean': 0,
            'seq_length_std': 0,
            'gc_content_mean': 0,
            'gc_content_std': 0,
            'mutation_load': 0  # Average mutations from consensus
        }
    
    # Basic diversity
    unique_seqs = len(set(sequences))
    seq_diversity = unique_seqs / len(sequences)
    
    # Similarity to consensus
    similarities = [calculate_sequence_similarity(seq, consensus) for seq in sequences]
    sim_mean = np.mean(similarities)
    sim_std = np.std(similarities) if len(similarities) > 1 else 0
    
    # Length statistics
    lengths = [len(seq) for seq in sequences]
    length_mean = np.mean(lengths)
    length_std = np.std(lengths) if len(lengths) > 1 else 0
    
    # GC content
    gc_contents = []
    for seq in sequences:
        seq_upper = seq.upper()
        gc_count = seq_upper.count('G') + seq_upper.count('C')
        gc_content = gc_count / len(seq) if seq else 0
        gc_contents.append(gc_content)
    
    gc_mean = np.mean(gc_contents)
    gc_std = np.std(gc_contents) if len(gc_contents) > 1 else 0
    
    # Mutation load (average number of mismatches from consensus)
    mutation_loads = []
    for seq in sequences:
        if len(seq) == len(consensus):
            mismatches = sum(1 for i, (a, b) in enumerate(zip(seq.upper(), consensus.upper())) 
                           if a != b and a in 'ATGC' and b in 'ATGC')
            mutation_loads.append(mismatches / len(seq))
    
    mutation_load = np.mean(mutation_loads) if mutation_loads else 0
    
    return {
        'seq_diversity': seq_diversity,
        'consensus_similarity_mean': sim_mean,
        'consensus_similarity_std': sim_std,
        'seq_length_mean': length_mean,
        'seq_length_std': length_std,
        'gc_content_mean': gc_mean,
        'gc_content_std': gc_std,
        'mutation_load': mutation_load
    }

def process_consolidated_motifs(tsv_file, expression_file=None, top_n_motifs=None, min_sites=10):
    """
    Process consolidated motif data into regression features
    
    Args:
        tsv_file: Path to consolidated_streme_sites.tsv
        expression_file: Optional path to expression data
        top_n_motifs: Only include top N most common motifs
        min_sites: Minimum sites required for a motif to be included
    
    Returns:
        Dictionary with feature matrix and metadata
    """
    
    print(f"Processing motif data from: {tsv_file}")
    
    # Read consolidated motif data
    line_gene_motifs = defaultdict(lambda: defaultdict(list))  # line -> gene -> list of motif data
    motif_global_stats = defaultdict(list)  # motif_id -> list of all sites
    
    try:
        with open(tsv_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                line = row['line']
                gene = row['gene_id']
                motif_id = row['consolidated_motif_id']
                position = int(row['start_pos'])
                sequence = row['sequence']
                strand = row['strand']
                
                motif_data = {
                    'motif_id': motif_id,
                    'position': position,
                    'sequence': sequence,
                    'strand': strand,
                    'consensus': row['motif_consensus']
                }
                
                line_gene_motifs[line][gene].append(motif_data)
                motif_global_stats[motif_id].append(motif_data)
    
    except Exception as e:
        print(f"Error reading motif file: {e}")
        sys.exit(1)
    
    # Filter motifs by frequency
    if top_n_motifs or min_sites:
        motif_counts = {motif: len(sites) for motif, sites in motif_global_stats.items()}
        
        if min_sites:
            motif_counts = {motif: count for motif, count in motif_counts.items() if count >= min_sites}
        
        if top_n_motifs:
            sorted_motifs = sorted(motif_counts.items(), key=lambda x: x[1], reverse=True)
            selected_motifs = set(motif for motif, count in sorted_motifs[:top_n_motifs])
        else:
            selected_motifs = set(motif_counts.keys())
        
        print(f"Selected {len(selected_motifs)} motifs for feature generation")
        print(f"Motif frequency range: {min(motif_counts.values())} - {max(motif_counts.values())} sites")
    
    # Generate features for each line-gene combination
    features = []
    feature_names = set()
    
    all_lines = sorted(line_gene_motifs.keys())
    all_genes = set()
    for line_data in line_gene_motifs.values():
        all_genes.update(line_data.keys())
    all_genes = sorted(all_genes)
    
    print(f"Generating features for {len(all_lines)} lines and {len(all_genes)} genes")
    print(f"This will create {len(all_lines) * len(all_genes)} gene×line samples")
    
    for line in all_lines:
        for gene in all_genes:
            gene_motifs = line_gene_motifs[line].get(gene, [])
            
            # Group motifs by type for this specific gene×line
            motifs_by_type = defaultdict(list)
            for motif_data in gene_motifs:
                motif_id = motif_data['motif_id']
                if not selected_motifs or motif_id in selected_motifs:
                    motifs_by_type[motif_id].append(motif_data)
            
            # Generate features for this gene×line combination
            row_features = {
                'line': line,
                'gene': gene,
                'gene_line_id': f"{gene}_{line}"  # Unique identifier
            }
            
            # For EVERY selected motif, generate features (even if not present)
            for motif_id in selected_motifs:
                motif_sites = motifs_by_type.get(motif_id, [])
                
                if motif_sites:
                    # Motif is present - calculate real features
                    positions = [site['position'] for site in motif_sites]
                    sequences = [site['sequence'] for site in motif_sites]
                    consensus = motif_sites[0]['consensus']
                    
                    # Calculate positional features
                    pos_features = calculate_position_features(positions)
                    for key, value in pos_features.items():
                        feature_name = f"{motif_id}_{key}"
                        row_features[feature_name] = value
                        feature_names.add(feature_name)
                    
                    # Calculate sequence variation features
                    var_features = calculate_sequence_variation_features(sequences, consensus)
                    for key, value in var_features.items():
                        feature_name = f"{motif_id}_{key}"
                        row_features[feature_name] = value
                        feature_names.add(feature_name)
                    
                    # Binary presence feature
                    feature_name = f"{motif_id}_present"
                    row_features[feature_name] = 1
                    feature_names.add(feature_name)
                
                else:
                    # Motif is NOT present - set all features to 0
                    for suffix in ['motif_count', 'position_mean', 'position_std', 'position_min', 
                                 'position_max', 'position_range', 'density_proximal', 'density_core', 
                                 'density_distal', 'spacing_regularity', 'seq_diversity', 
                                 'consensus_similarity_mean', 'consensus_similarity_std', 
                                 'seq_length_mean', 'seq_length_std', 'gc_content_mean', 
                                 'gc_content_std', 'mutation_load', 'present']:
                        feature_name = f"{motif_id}_{suffix}"
                        row_features[feature_name] = 0
                        feature_names.add(feature_name)
            
            features.append(row_features)
    
    # Load expression data if provided
    expression_data = {}
    if expression_file and os.path.exists(expression_file):
        print(f"Loading expression data from: {expression_file}")
        try:
            with open(expression_file, 'r') as f:
                reader = csv.DictReader(f, delimiter='\t')
                for row in reader:
                    line = row.get('line', row.get('Line', ''))
                    gene = row.get('gene', row.get('Gene', row.get('gene_id', '')))
                    expression = float(row.get('expression', row.get('Expression', row.get('value', 0))))
                    expression_data[(line, gene)] = expression
        except Exception as e:
            print(f"Warning: Could not load expression data: {e}")
    
    # Add expression values to features
    if expression_data:
        for row in features:
            key = (row['line'], row['gene'])
            row['expression'] = expression_data.get(key, None)
    
    return {
        'features': features,
        'feature_names': sorted(feature_names),
        'motif_stats': motif_global_stats,
        'selected_motifs': selected_motifs if 'selected_motifs' in locals() else set(motif_global_stats.keys()),
        'n_lines': len(all_lines),
        'n_genes': len(all_genes),
        'n_features': len(feature_names)
    }

def write_features_to_files(results, output_prefix):
    """Write feature matrix to CSV files"""
    
    features = results['features']
    feature_names = results['feature_names']
    
    # Write main feature matrix
    feature_file = f"{output_prefix}_features.csv"
    
    # Ensure consistent column ordering: identifiers first, then all motif features, then expression
    all_columns = ['line', 'gene', 'gene_line_id'] + sorted(feature_names)
    if any('expression' in row for row in features):
        all_columns.append('expression')
    
    with open(feature_file, 'w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=all_columns)
        writer.writeheader()
        for row in features:
            # Ensure all feature columns are present (fill missing with 0)
            complete_row = {col: row.get(col, 0) for col in all_columns}
            # But keep string columns as strings
            for id_col in ['line', 'gene', 'gene_line_id']:
                if id_col in row:
                    complete_row[id_col] = row[id_col]
            writer.writerow(complete_row)
    
    print(f"Feature matrix written to: {feature_file}")
    print(f"Shape: {len(features)} samples × {len(feature_names)} features")
    print(f"Each row = one gene×line combination with features for ALL motifs")
    
    # Write feature descriptions
    desc_file = f"{output_prefix}_feature_descriptions.txt"
    with open(desc_file, 'w') as f:
        f.write("# Motif-based Regression Features\n\n")
        f.write(f"Generated from motif data with {results['n_lines']} lines, {results['n_genes']} genes\n")
        f.write(f"Selected {len(results['selected_motifs'])} motifs for feature generation\n\n")
        
        f.write("## Data Format:\n")
        f.write("- Each row represents one gene×line combination\n")
        f.write("- Each gene appears once per line (even if no motifs present)\n")
        f.write("- Features for ALL selected motifs are included (0 if motif absent)\n")
        f.write("- Ready for regression with expression as target variable\n\n")
        
        f.write("## Columns:\n")
        f.write("- line: Line identifier\n")
        f.write("- gene: Gene identifier\n") 
        f.write("- gene_line_id: Unique identifier (gene_line)\n")
        f.write("- [motif_features]: Features for each motif (see categories below)\n")
        f.write("- expression: Target variable (if provided)\n\n")
        
        f.write("## Feature Categories:\n\n")
        f.write("### Positional Features (per motif):\n")
        f.write("- motif_count: Number of motif occurrences\n")
        f.write("- position_mean/std/min/max/range: Position statistics\n")
        f.write("- density_proximal: Motif density close to TSS (700-1000bp upstream)\n")
        f.write("- density_core: Motif density in middle region (300-700bp upstream)\n")
        f.write("- density_distal: Motif density far from TSS (0-300bp upstream)\n")
        f.write("- spacing_regularity: Regularity of motif spacing (lower = more regular)\n\n")
        
        f.write("### Sequence Variation Features (per motif):\n")
        f.write("- seq_diversity: Fraction of unique sequences\n")
        f.write("- consensus_similarity_mean/std: Similarity to consensus sequence\n")
        f.write("- seq_length_mean/std: Length variation\n")
        f.write("- gc_content_mean/std: GC content variation\n")
        f.write("- mutation_load: Average mutations from consensus\n\n")
        
        f.write("### Presence Features (per motif):\n")
        f.write("- present: Binary indicator of motif presence\n\n")
        
        f.write("## Selected Motifs:\n")
        for motif in sorted(results['selected_motifs']):
            count = len(results['motif_stats'][motif])
            f.write(f"- {motif}: {count} total sites\n")
    
    print(f"Feature descriptions written to: {desc_file}")
    
    # Write summary statistics
    summary_file = f"{output_prefix}_summary.txt"
    with open(summary_file, 'w') as f:
        f.write("# Regression Feature Generation Summary\n\n")
        f.write(f"Total samples (gene×line combinations): {len(features)}\n")
        f.write(f"Total features: {len(feature_names)}\n")
        f.write(f"Lines: {results['n_lines']}\n")
        f.write(f"Genes: {results['n_genes']}\n")
        f.write(f"Motifs included: {len(results['selected_motifs'])}\n\n")
        
        f.write("## Matrix Structure:\n")
        f.write("- Rows: Each gene×line combination (one row per gene per line)\n")
        f.write("- Columns: Features for ALL motifs (present=real values, absent=0)\n")
        f.write("- Format: Ready for ML regression (X matrix for predicting expression)\n\n")
        
        # Feature type breakdown
        pos_features = sum(1 for name in feature_names if any(suffix in name for suffix in 
                          ['count', 'position', 'density', 'spacing']))
        var_features = sum(1 for name in feature_names if any(suffix in name for suffix in 
                          ['diversity', 'similarity', 'length', 'gc_content', 'mutation']))
        pres_features = sum(1 for name in feature_names if 'present' in name)
        
        f.write(f"Feature breakdown:\n")
        f.write(f"- Positional features: {pos_features}\n")
        f.write(f"- Sequence variation features: {var_features}\n")
        f.write(f"- Presence features: {pres_features}\n\n")
        
        f.write("## Usage Notes:\n")
        f.write("1. Features are designed for regression/ML models predicting gene expression\n")
        f.write("2. Consider feature scaling/normalization before modeling\n")
        f.write("3. High-dimensional data - consider feature selection or regularization\n")
        f.write("4. Check for multicollinearity between related features\n")
        f.write("5. Positional features assume 1000bp upstream promoter sequences\n")
        f.write("6. Position 0 = far from TSS, Position 1000 = close to TSS\n")
    
    print(f"Summary written to: {summary_file}")

def main():
    parser = argparse.ArgumentParser(
        description="Convert consolidated motif data into regression features",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic feature generation
  python motif_to_regression_features.py consolidated_streme_sites.tsv

  # Include expression data and select top 100 motifs
  python motif_to_regression_features.py consolidated_streme_sites.tsv \\
    --expression expression_data.tsv --top-motifs 100

  # Filter motifs with minimum 50 sites
  python motif_to_regression_features.py consolidated_streme_sites.tsv \\
    --min-sites 50 --output-prefix my_analysis
        """
    )
    
    parser.add_argument('motif_file', help='Consolidated motif TSV file')
    parser.add_argument('--expression', '-e', help='Expression data file (TSV with line, gene, expression columns)')
    parser.add_argument('--top-motifs', '-t', type=int, help='Include only top N most frequent motifs')
    parser.add_argument('--min-sites', '-m', type=int, default=10, help='Minimum sites required for motif inclusion (default: 10)')
    parser.add_argument('--output-prefix', '-o', default='motif_regression', help='Output file prefix (default: motif_regression)')
    
    args = parser.parse_args()
    
    # Process motif data
    results = process_consolidated_motifs(
        args.motif_file,
        expression_file=args.expression,
        top_n_motifs=args.top_motifs,
        min_sites=args.min_sites
    )
    
    # Write output files
    write_features_to_files(results, args.output_prefix)
    
    print(f"\n✅ Feature generation complete!")
    print(f"Ready for regression modeling with {results['n_features']} features")

if __name__ == "__main__":
    main()
