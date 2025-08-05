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
import csv
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

def extract_line_name(directory_name):
    """
    Extract clean line name from STREME directory name
    Examples:
    - streme_Genes_IM155_DNA_lifted -> IM155
    - streme_IM502 -> IM502
    """
    # Remove 'streme_' prefix
    name = directory_name.replace('streme_', '')
    
    # Extract IM number using regex
    match = re.search(r'IM\d+', name)
    if match:
        return match.group()
    
    # Fallback: return the cleaned name
    return name

def parse_sites_file(sites_file, line_name):
    """
    Parse a STREME sites.tsv file and return structured data
    """
    print(f"  Processing {sites_file} for line {line_name}")
    
    try:
        # Read the TSV file manually to handle NaN values better
        sites_data = []
        
        with open(sites_file, 'r') as f:
            lines = f.readlines()
        
        if not lines:
            print(f"    Warning: Empty file {sites_file}")
            return []
        
        # Parse header
        header = lines[0].strip().split('\t')
        
        # Validate required columns
        required_cols = ['motif_ID', 'seq_ID', 'site_Start', 'site_End', 'site_Strand', 'site_Score', 'site_Sequence']
        col_indices = {}
        
        for col in required_cols:
            if col in header:
                col_indices[col] = header.index(col)
            else:
                print(f"Warning: Missing column '{col}' in {sites_file}")
                return []
        
        # Parse data rows
        for line_num, line in enumerate(lines[1:], 2):
            if not line.strip():
                continue
            
            # Skip comment lines (STREME motif separators)
            if line.strip().startswith('#'):
                continue
                
            fields = line.strip().split('\t')
            
            try:
                # Check if we have enough fields
                if len(fields) < len(col_indices):
                    print(f"    Skipping line {line_num}: insufficient columns ({len(fields)} < {len(col_indices)})")
                    continue
                
                # Check for missing or NaN values
                start_val = fields[col_indices['site_Start']]
                end_val = fields[col_indices['site_End']]
                score_val = fields[col_indices['site_Score']]
                sequence_val = fields[col_indices['site_Sequence']]
                motif_id_val = fields[col_indices['motif_ID']]
                
                # Skip rows with missing essential data
                if (start_val in ['', 'NaN', 'nan'] or 
                    end_val in ['', 'NaN', 'nan'] or 
                    score_val in ['', 'NaN', 'nan'] or 
                    sequence_val in ['', 'NaN', 'nan'] or
                    motif_id_val in ['', 'NaN', 'nan']):
                    print(f"    Skipping line {line_num}: missing essential data")
                    continue
                
                site_data = {
                    'original_motif_id': motif_id_val,
                    'line': line_name,
                    'gene_id': fields[col_indices['seq_ID']],
                    'start_pos': int(float(start_val)),  # Convert via float first to handle decimal strings
                    'end_pos': int(float(end_val)),
                    'strand': fields[col_indices['site_Strand']],
                    'score': float(score_val),
                    'sequence': sequence_val.upper(),
                    'motif_alt_id': fields[col_indices.get('motif_ALT_ID', 0)] if 'motif_ALT_ID' in col_indices else '',
                    'length': len(sequence_val)
                }
                sites_data.append(site_data)
                
            except (ValueError, IndexError) as e:
                print(f"    Skipping line {line_num}: {e}")
                continue
        
        print(f"    Found {len(sites_data)} valid motif sites")
        return sites_data
        
    except Exception as e:
        print(f"Error processing {sites_file}: {e}")
        return []

def extract_motif_consensus(motif_id):
    """
    Extract the consensus sequence from motif_ID
    Examples:
    - 1-AAAARAAAAAAAAAA -> AAAARAAAAAAAAAA
    - 2-CCCCGTTTTCCCCCC -> CCCCGTTTTCCCCCC
    """
    # Remove numerical prefix (e.g., "1-", "2-", etc.)
    if '-' in motif_id:
        return motif_id.split('-', 1)[1]
    return motif_id

def consolidate_motifs_by_sequence(all_sites, similarity_threshold=0.75):
    """
    Consolidate motifs based on motif_ID consensus patterns, not actual sequences found
    """
    print(f"Consolidating motifs using similarity threshold {similarity_threshold}")
    
    # Collect all unique motif consensus patterns
    unique_motifs = {}
    
    for site in all_sites:
        motif_consensus = extract_motif_consensus(site['original_motif_id'])
        if motif_consensus not in unique_motifs:
            unique_motifs[motif_consensus] = []
        unique_motifs[motif_consensus].append(site)
    
    print(f"Found {len(unique_motifs)} unique motif patterns before consolidation")
    
    # Cluster similar motif consensus patterns
    clusters = []
    processed_motifs = set()
    
    for motif1 in unique_motifs:
        if motif1 in processed_motifs:
            continue
            
        # Start new cluster
        cluster_motifs = [motif1]
        processed_motifs.add(motif1)
        
        # Find similar motif patterns
        for motif2 in unique_motifs:
            if motif2 in processed_motifs:
                continue
                
            if sequences_are_similar(motif1, motif2, similarity_threshold):
                cluster_motifs.append(motif2)
                processed_motifs.add(motif2)
        
        clusters.append(cluster_motifs)
    
    print(f"Consolidated into {len(clusters)} motif clusters")
    
    # Assign consolidated IDs
    consolidated_data = []
    
    for cluster_idx, cluster_motifs in enumerate(clusters):
        consolidated_id = f"MOTIF_{cluster_idx + 1:03d}"
        
        # Get representative motif (most common one in cluster)
        motif_counts = {}
        for motif in cluster_motifs:
            motif_counts[motif] = len(unique_motifs[motif])
        
        representative_motif = max(motif_counts, key=motif_counts.get)
        
        # Add all sites from this cluster
        cluster_sites = []
        for motif in cluster_motifs:
            for site in unique_motifs[motif]:
                site_copy = site.copy()
                site_copy['consolidated_motif_id'] = consolidated_id
                site_copy['representative_sequence'] = representative_motif
                site_copy['cluster_size'] = len(cluster_motifs)
                site_copy['motif_consensus'] = motif  # Add the consensus pattern
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
        'consolidated_motif_id', 'original_motif_id', 'motif_consensus', 'line', 'gene_id',
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
            # Extract line name properly
            line_name = extract_line_name(item.name)
            
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
