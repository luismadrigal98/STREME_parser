#!/usr/bin/env python3
"""
STREME Sites Consolidator - Parse sites.tsv files and create comprehensive motif maps

This tool processes STREME sites.tsv files from multiple lines, consolidates
similar motifs using IUPAC-aware similarity, and creates comprehensive
gene-level regulatory maps with exact positions and sequences.
"""

import os
import re
import sys
import argparse
import pandas as pd
from pathlib import Path
from collections import defaultdict

def expand_iupac_for_comparison(seq1, seq2):
    """
    Expand IUPAC codes and find best possible alignment score.
    Returns normalized similarity score (0-1).
    """
    iupac_map = {
        'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C',
        'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT',
        'K': 'GT', 'M': 'AC', 'B': 'CGT', 'D': 'AGT',
        'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'
    }
    
    def get_possible_bases(char):
        return iupac_map.get(char.upper(), char.upper())
    
    # Calculate maximum possible matches by position
    max_matches = 0
    total_positions = min(len(seq1), len(seq2))
    
    for i in range(total_positions):
        bases1 = set(get_possible_bases(seq1[i]))
        bases2 = set(get_possible_bases(seq2[i]))
        if bases1 & bases2:  # Any overlap
            max_matches += 1
    
    # Penalize length differences heavily
    len_diff = abs(len(seq1) - len(seq2))
    len_penalty = len_diff * 0.2  # Heavy penalty for length differences
    
    if total_positions == 0:
        return 0.0
    
    # Calculate similarity and apply length penalty
    base_similarity = max_matches / total_positions
    final_similarity = base_similarity - len_penalty
    
    return max(0.0, final_similarity)

def sequences_are_similar(seq1, seq2, threshold=0.75):
    """
    Check if two sequences are similar enough to be the same motif.
    Uses proper sequence comparison with IUPAC handling and length penalties.
    """
    # Quick filters
    if abs(len(seq1) - len(seq2)) > 3:  # Too different in length
        return False
    
    if len(seq1) < 3 or len(seq2) < 3:  # Too short to compare reliably
        return False
    
    # Use IUPAC-aware comparison with length penalties
    similarity = expand_iupac_for_comparison(seq1, seq2)
    
    return similarity >= threshold

def parse_sites_file(sites_file, line_name):
    """
    Parse a STREME sites.tsv file and return structured data
    """
    print(f"  Processing {sites_file} for line {line_name}")
    
    try:
        # Read the TSV file
        df = pd.read_csv(sites_file, sep='\t')
        
        # Validate required columns
        required_cols = ['motif_ID', 'seq_ID', 'site_Start', 'site_End', 'site_Strand', 'site_Score', 'site_Sequence']
        missing_cols = [col for col in required_cols if col not in df.columns]
        if missing_cols:
            print(f"Warning: Missing columns in {sites_file}: {missing_cols}")
            return []
        
        sites_data = []
        
        for _, row in df.iterrows():
            site_data = {
                'original_motif_id': row['motif_ID'],
                'line': line_name,
                'gene_id': row['seq_ID'],
                'start_pos': int(row['site_Start']),
                'end_pos': int(row['site_End']),
                'strand': row['site_Strand'],
                'score': float(row['site_Score']),
                'sequence': row['site_Sequence'].upper(),
                'motif_alt_id': row.get('motif_ALT_ID', ''),
                'length': len(row['site_Sequence'])
            }
            sites_data.append(site_data)
        
        print(f"    Found {len(sites_data)} motif sites")
        return sites_data
        
    except Exception as e:
        print(f"Error processing {sites_file}: {e}")
        return []

def consolidate_motifs_by_sequence(all_sites, similarity_threshold=0.75):
    """
    Consolidate motifs based on actual sequences found, not just consensus
    """
    print(f"Consolidating motifs using similarity threshold {similarity_threshold}")
    
    # Collect all unique sequences
    unique_sequences = {}
    sequence_to_cluster = {}
    
    for site in all_sites:
        seq = site['sequence']
        if seq not in unique_sequences:
            unique_sequences[seq] = []
        unique_sequences[seq].append(site)
    
    print(f"Found {len(unique_sequences)} unique sequences before consolidation")
    
    # Cluster similar sequences
    clusters = []
    processed_sequences = set()
    
    for seq1 in unique_sequences:
        if seq1 in processed_sequences:
            continue
            
        # Start new cluster
        cluster_sequences = [seq1]
        processed_sequences.add(seq1)
        
        # Find similar sequences
        for seq2 in unique_sequences:
            if seq2 in processed_sequences:
                continue
                
            if sequences_are_similar(seq1, seq2, similarity_threshold):
                cluster_sequences.append(seq2)
                processed_sequences.add(seq2)
        
        clusters.append(cluster_sequences)
    
    print(f"Consolidated into {len(clusters)} motif clusters")
    
    # Assign consolidated IDs
    consolidated_data = []
    
    for cluster_idx, cluster_seqs in enumerate(clusters):
        consolidated_id = f"MOTIF_{cluster_idx + 1:03d}"
        
        # Get representative sequence (most common one in cluster)
        seq_counts = {}
        for seq in cluster_seqs:
            seq_counts[seq] = len(unique_sequences[seq])
        
        representative_seq = max(seq_counts, key=seq_counts.get)
        
        # Add all sites from this cluster
        cluster_sites = []
        for seq in cluster_seqs:
            for site in unique_sequences[seq]:
                site_copy = site.copy()
                site_copy['consolidated_motif_id'] = consolidated_id
                site_copy['representative_sequence'] = representative_seq
                site_copy['cluster_size'] = len(cluster_seqs)
                cluster_sites.append(site_copy)
        
        consolidated_data.extend(cluster_sites)
    
    return consolidated_data

def add_relative_positions(consolidated_data):
    """
    Add relative position information per gene per line
    """
    print("Adding relative position information...")
    
    # Group by gene and line
    gene_line_groups = defaultdict(list)
    
    for site in consolidated_data:
        key = (site['gene_id'], site['line'])
        gene_line_groups[key].append(site)
    
    # Sort by position and add relative position info
    for (gene_id, line), sites in gene_line_groups.items():
        # Sort by start position
        sites.sort(key=lambda x: x['start_pos'])
        
        # Add relative position information
        for idx, site in enumerate(sites):
            site['relative_position'] = idx + 1  # 1-based indexing
            site['total_motifs_in_gene'] = len(sites)
            site['relative_position_fraction'] = (idx + 1) / len(sites)
    
    return consolidated_data

def write_comprehensive_output(consolidated_data, output_file):
    """
    Write comprehensive motif mapping results
    """
    print(f"Writing comprehensive results to {output_file}")
    
    # Convert to DataFrame for easy manipulation
    df = pd.DataFrame(consolidated_data)
    
    # Sort by consolidated motif ID, then by line, then by gene, then by position
    df = df.sort_values(['consolidated_motif_id', 'line', 'gene_id', 'start_pos'])
    
    # Select and order columns
    output_columns = [
        'consolidated_motif_id', 'original_motif_id', 'line', 'gene_id',
        'start_pos', 'end_pos', 'strand', 'score', 'sequence',
        'representative_sequence', 'relative_position', 'total_motifs_in_gene',
        'relative_position_fraction', 'cluster_size', 'length'
    ]
    
    # Ensure all columns exist
    for col in output_columns:
        if col not in df.columns:
            df[col] = ''
    
    df[output_columns].to_csv(output_file, sep='\t', index=False)
    
    return df

def write_summary_statistics(df, summary_file):
    """
    Write summary statistics about the consolidation
    """
    print(f"Writing summary statistics to {summary_file}")
    
    with open(summary_file, 'w') as f:
        f.write("# STREME Sites Consolidation Summary\n\n")
        
        # Basic statistics
        f.write(f"Total motif sites: {len(df)}\n")
        f.write(f"Unique consolidated motifs: {df['consolidated_motif_id'].nunique()}\n")
        f.write(f"Lines analyzed: {df['line'].nunique()}\n")
        f.write(f"Genes with motifs: {df['gene_id'].nunique()}\n\n")
        
        # Per-line statistics
        f.write("## Statistics by Line\n")
        line_stats = df.groupby('line').agg({
            'gene_id': 'nunique',
            'consolidated_motif_id': 'nunique',
            'sequence': 'count'
        }).rename(columns={
            'gene_id': 'Genes_with_motifs',
            'consolidated_motif_id': 'Unique_motifs',
            'sequence': 'Total_sites'
        })
        f.write(line_stats.to_string())
        f.write("\n\n")
        
        # Motif frequency
        f.write("## Top 20 Most Common Motifs\n")
        motif_freq = df['consolidated_motif_id'].value_counts().head(20)
        for motif, count in motif_freq.items():
            repr_seq = df[df['consolidated_motif_id'] == motif]['representative_sequence'].iloc[0]
            f.write(f"{motif}\t{repr_seq}\t{count} sites\n")

def main():
    parser = argparse.ArgumentParser(
        description='Consolidate STREME sites.tsv files across multiple lines',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic consolidation
  %(prog)s /path/to/streme/results/ --output consolidated_sites
  
  # With custom similarity threshold
  %(prog)s /path/to/streme/results/ --threshold 0.8 --output results/
  
  # Process specific lines only
  %(prog)s /path/to/streme/results/ --lines IM502,IM664 --output targeted/

Input Directory Structure:
  streme_results/
  ├── streme_IM502/sites.tsv
  ├── streme_IM664/sites.tsv
  └── streme_IM767/sites.tsv
        """
    )
    
    parser.add_argument(
        'streme_results_dir',
        help='Directory containing STREME output folders with sites.tsv files'
    )
    
    parser.add_argument(
        '--output', '-o',
        default='consolidated_streme_sites',
        help='Output file prefix (default: consolidated_streme_sites)'
    )
    
    parser.add_argument(
        '--threshold', '-t',
        type=float,
        default=0.75,
        help='Similarity threshold for motif consolidation (default: 0.75)'
    )
    
    parser.add_argument(
        '--lines', '-l',
        help='Comma-separated list of specific lines to process (e.g., IM502,IM664)'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output'
    )
    
    args = parser.parse_args()
    
    # Validate input directory
    results_dir = Path(args.streme_results_dir)
    if not results_dir.exists():
        print(f"Error: Directory {results_dir} does not exist")
        sys.exit(1)
    
    print("=== STREME Sites Consolidator ===")
    print(f"Input directory: {results_dir}")
    print(f"Similarity threshold: {args.threshold}")
    print(f"Output prefix: {args.output}")
    
    # Find STREME directories with sites.tsv files
    streme_dirs = []
    target_lines = args.lines.split(',') if args.lines else None
    
    for item in results_dir.iterdir():
        if item.is_dir() and item.name.startswith('streme_'):
            # Extract line name
            line_name = item.name.replace('streme_', '')
            
            # Skip if not in target lines
            if target_lines and line_name not in target_lines:
                continue
            
            sites_file = item / 'sites.tsv'
            if sites_file.exists():
                streme_dirs.append((sites_file, line_name))
            elif args.verbose:
                print(f"Warning: No sites.tsv found in {item}")
    
    if not streme_dirs:
        print("Error: No STREME directories with sites.tsv files found")
        sys.exit(1)
    
    print(f"Found {len(streme_dirs)} STREME directories to process")
    
    # Parse all sites files
    all_sites = []
    for sites_file, line_name in streme_dirs:
        sites_data = parse_sites_file(sites_file, line_name)
        all_sites.extend(sites_data)
    
    if not all_sites:
        print("Error: No motif sites found in any file")
        sys.exit(1)
    
    print(f"Total sites collected: {len(all_sites)}")
    
    # Consolidate motifs
    consolidated_data = consolidate_motifs_by_sequence(all_sites, args.threshold)
    
    # Add relative position information
    consolidated_data = add_relative_positions(consolidated_data)
    
    # Write outputs
    output_file = f"{args.output}.tsv"
    summary_file = f"{args.output}_summary.txt"
    
    df = write_comprehensive_output(consolidated_data, output_file)
    write_summary_statistics(df, summary_file)
    
    print(f"\n=== COMPLETED ===")
    print(f"Comprehensive results: {output_file}")
    print(f"Summary statistics: {summary_file}")
    print(f"Total consolidated motifs: {df['consolidated_motif_id'].nunique()}")
    print(f"Total sites: {len(df)}")

if __name__ == "__main__":
    main()
