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
from pathlib import Path
from collections import defaultdict

try:
    import pandas as pd
    HAS_PANDAS = True
except ImportError:
    HAS_PANDAS = False

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

def calculate_consensus_from_sequences(sequences):
    """
    Calculate a true consensus sequence from actual sequences found,
    not just using the STREME consensus pattern.
    """
    if not sequences:
        return ""
    
    # Find the most common length
    lengths = [len(seq) for seq in sequences]
    most_common_length = max(set(lengths), key=lengths.count)
    
    # Filter to sequences of the most common length
    filtered_seqs = [seq for seq in sequences if len(seq) == most_common_length]
    
    if not filtered_seqs:
        return sequences[0]  # Fallback
    
    # Calculate position-wise consensus
    consensus = []
    for pos in range(most_common_length):
        bases_at_pos = [seq[pos] for seq in filtered_seqs if pos < len(seq)]
        
        # Count each base
        base_counts = {'A': 0, 'T': 0, 'G': 0, 'C': 0}
        for base in bases_at_pos:
            if base in base_counts:
                base_counts[base] += 1
        
        # Find most common base or use IUPAC code
        total_bases = sum(base_counts.values())
        if total_bases == 0:
            consensus.append('N')
            continue
            
        # Calculate percentages
        percentages = {base: count/total_bases for base, count in base_counts.items()}
        
        # Use IUPAC codes for ambiguous positions
        if percentages['A'] >= 0.8:
            consensus.append('A')
        elif percentages['T'] >= 0.8:
            consensus.append('T')
        elif percentages['G'] >= 0.8:
            consensus.append('G')
        elif percentages['C'] >= 0.8:
            consensus.append('C')
        elif percentages['A'] + percentages['T'] >= 0.8:
            consensus.append('W')  # A or T
        elif percentages['G'] + percentages['C'] >= 0.8:
            consensus.append('S')  # G or C
        elif percentages['A'] + percentages['G'] >= 0.8:
            consensus.append('R')  # A or G
        elif percentages['C'] + percentages['T'] >= 0.8:
            consensus.append('Y')  # C or T
        else:
            consensus.append('N')  # Highly variable
    
    return ''.join(consensus)

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
    Three-pass consolidation based on TRUE consensus patterns calculated from actual sequences:
    1. Calculate true consensus for each STREME motif pattern
    2. Cluster STREME motifs based on true consensus similarity  
    3. Calculate final consensus from all sequences in each cluster
    """
    print(f"Three-pass consolidation using similarity threshold {similarity_threshold}")
    
    # PASS 1: Calculate true consensus for each original STREME motif
    print("Pass 1: Calculating true consensus for each STREME motif...")
    unique_motifs = {}
    motif_true_consensus = {}
    
    for site in all_sites:
        motif_consensus = extract_motif_consensus(site['original_motif_id'])
        if motif_consensus not in unique_motifs:
            unique_motifs[motif_consensus] = []
        unique_motifs[motif_consensus].append(site)
    
    # Calculate true consensus for each STREME motif
    for motif_pattern, sites in unique_motifs.items():
        sequences = [site['sequence'] for site in sites]
        true_consensus = calculate_consensus_from_sequences(sequences)
        motif_true_consensus[motif_pattern] = true_consensus
        print(f"  {motif_pattern} -> {true_consensus} (from {len(sequences)} sequences)")
    
    print(f"Found {len(unique_motifs)} unique STREME motif patterns")
    
    # PASS 2: Cluster based on TRUE consensus patterns
    print("Pass 2: Clustering based on true consensus patterns...")
    clusters = []
    processed_motifs = set()
    
    for motif1 in unique_motifs:
        if motif1 in processed_motifs:
            continue
            
        # Start new cluster
        cluster_motifs = [motif1]
        processed_motifs.add(motif1)
        true_consensus_1 = motif_true_consensus[motif1]
        
        # Find similar TRUE consensus patterns
        for motif2 in unique_motifs:
            if motif2 in processed_motifs:
                continue
            
            true_consensus_2 = motif_true_consensus[motif2]
            
            if sequences_are_similar(true_consensus_1, true_consensus_2, similarity_threshold):
                cluster_motifs.append(motif2)
                processed_motifs.add(motif2)
                print(f"  Clustered {motif1} ({true_consensus_1}) with {motif2} ({true_consensus_2})")
        
        clusters.append(cluster_motifs)
    
    print(f"Consolidated into {len(clusters)} motif clusters based on true consensus")
    
    # PASS 3: Assign consolidated IDs and calculate final true consensus per cluster
    consolidated_data = []
    
    for cluster_idx, cluster_motifs in enumerate(clusters):
        consolidated_id = f"MOTIF_{cluster_idx + 1:03d}"
        
        # Collect all actual sequences from this cluster
        all_cluster_sequences = []
        cluster_sites = []
        
        for motif in cluster_motifs:
            for site in unique_motifs[motif]:
                all_cluster_sequences.append(site['sequence'])
                site_copy = site.copy()
                site_copy['consolidated_motif_id'] = consolidated_id
                site_copy['cluster_size'] = len(cluster_motifs)
                site_copy['original_streme_consensus'] = motif  # Keep original STREME consensus
                cluster_sites.append(site_copy)
        
        # Calculate final true consensus from ALL sequences in the cluster
        final_true_consensus = calculate_consensus_from_sequences(all_cluster_sequences)
        
        # Get most common original motif pattern for reference
        motif_counts = {}
        for motif in cluster_motifs:
            motif_counts[motif] = len(unique_motifs[motif])
        most_common_original = max(motif_counts, key=motif_counts.get)
        
        # Update all sites in this cluster
        for site in cluster_sites:
            site['motif_consensus'] = final_true_consensus  # True consensus from actual sequences
            site['original_streme_consensus'] = motif  # Original STREME pattern for reference
        
        if len(cluster_motifs) > 1:
            print(f"  Cluster {consolidated_id}: {len(cluster_motifs)} STREME motifs -> {len(all_cluster_sequences)} sequences -> final consensus: {final_true_consensus}")
        
        consolidated_data.extend(cluster_sites)
    
    return consolidated_data

def merge_overlapping_motifs(motif_sites, overlap_threshold=0.5):
    """
    Merge overlapping motif occurrences of the same type within the same gene/line.
    This prevents STREME's sliding window artifacts from creating redundant entries.
    
    Args:
        motif_sites: List of motif site dictionaries
        overlap_threshold: Minimum fraction of overlap to consider for merging (0.5 = 50%)
    
    Returns:
        List of merged motif sites
    """
    print(f"Merging overlapping motifs with threshold {overlap_threshold}")
    
    # Group by gene, line, strand, and consolidated motif ID
    groups = {}
    
    for site in motif_sites:
        key = (site['gene_id'], site['line'], site['strand'], site['consolidated_motif_id'])
        if key not in groups:
            groups[key] = []
        groups[key].append(site)
    
    merged_sites = []
    total_original = len(motif_sites)
    total_merged = 0
    
    for key, sites in groups.items():
        # Sort by start position
        sites.sort(key=lambda x: x['start_pos'])
        
        # Merge overlapping sites
        merged_group = []
        
        for site in sites:
            if not merged_group:
                # First site in group
                site['merged_count'] = 1
                merged_group.append(site)
            else:
                last_merged = merged_group[-1]
                
                # Check for overlap
                overlap_start = max(last_merged['start_pos'], site['start_pos'])
                overlap_end = min(last_merged['end_pos'], site['end_pos'])
                overlap_length = max(0, overlap_end - overlap_start + 1)
                
                # Calculate overlap fraction relative to shorter sequence
                min_length = min(
                    last_merged['end_pos'] - last_merged['start_pos'] + 1,
                    site['end_pos'] - site['start_pos'] + 1
                )
                overlap_fraction = overlap_length / min_length if min_length > 0 else 0
                
                if overlap_fraction >= overlap_threshold:
                    # Merge with the last site
                    total_merged += 1
                    
                    # Extend the boundaries
                    last_merged['start_pos'] = min(last_merged['start_pos'], site['start_pos'])
                    last_merged['end_pos'] = max(last_merged['end_pos'], site['end_pos'])
                    
                    # Keep the best score
                    last_merged['score'] = max(last_merged['score'], site['score'])
                    
                    # Update sequence to cover the full merged region
                    # Keep the sequence with the better score
                    if site['score'] > last_merged['score'] or len(site['sequence']) > len(last_merged['sequence']):
                        last_merged['sequence'] = site['sequence']
                    
                    # Update length
                    last_merged['length'] = last_merged['end_pos'] - last_merged['start_pos'] + 1
                    
                    # Track merged count
                    last_merged['merged_count'] = last_merged.get('merged_count', 1) + 1
                    
                else:
                    # No overlap, add as separate site
                    site['merged_count'] = 1
                    merged_group.append(site)
        
        # Add merged sites to final list
        merged_sites.extend(merged_group)
    
    print(f"  Merged {total_merged} overlapping sites: {total_original} -> {len(merged_sites)} final sites")
    return merged_sites

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
    
    if HAS_PANDAS:
        # Use pandas for sorting and manipulation if available
        df = pd.DataFrame(consolidated_data)
        
        # Sort by consolidated motif ID, then by line, then by gene, then by position
        df = df.sort_values(['consolidated_motif_id', 'line', 'gene_id', 'start_pos'])
        
        # Select and order columns
        output_columns = [
            'consolidated_motif_id', 'original_motif_id', 'motif_consensus', 'original_streme_consensus', 
            'line', 'gene_id', 'start_pos', 'end_pos', 'strand', 'score', 'sequence',
            'relative_position', 'total_motifs_in_gene', 'relative_position_fraction', 
            'cluster_size', 'merged_count', 'length'
        ]
        
        # Ensure all columns exist
        for col in output_columns:
            if col not in df.columns:
                df[col] = ''
        
        df[output_columns].to_csv(output_file, sep='\t', index=False)
        return df
    else:
        # Fallback to manual CSV writing without pandas
        # Sort data manually
        consolidated_data.sort(key=lambda x: (x.get('consolidated_motif_id', ''), 
                                            x.get('line', ''), 
                                            x.get('gene_id', ''), 
                                            x.get('start_pos', 0)))
        
        output_columns = [
            'consolidated_motif_id', 'original_motif_id', 'motif_consensus', 'original_streme_consensus', 
            'line', 'gene_id', 'start_pos', 'end_pos', 'strand', 'score', 'sequence',
            'relative_position', 'total_motifs_in_gene', 'relative_position_fraction', 
            'cluster_size', 'merged_count', 'length'
        ]
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=output_columns, delimiter='\t')
            writer.writeheader()
            for row in consolidated_data:
                # Ensure all columns exist in row
                output_row = {col: row.get(col, '') for col in output_columns}
                writer.writerow(output_row)
        
        # Return simple dict-like object for consistency
        class SimpleDataFrame:
            def __init__(self, data):
                self.data = data
            
            def nunique(self, column):
                values = set(row.get(column) for row in self.data if row.get(column))
                return len(values)
            
            def __len__(self):
                return len(self.data)
        
        return SimpleDataFrame(consolidated_data)

def write_summary_statistics(df, summary_file):
    """
    Write summary statistics about the consolidation
    """
    print(f"Writing summary statistics to {summary_file}")
    
    with open(summary_file, 'w') as f:
        f.write("# STREME Sites Consolidation Summary\n\n")
        
        if HAS_PANDAS and hasattr(df, 'groupby'):
            # Use pandas functionality
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
                true_consensus = df[df['consolidated_motif_id'] == motif]['motif_consensus'].iloc[0]
                f.write(f"{motif}\t{true_consensus}\t{count} sites\n")
        else:
            # Manual statistics calculation
            data = df.data if hasattr(df, 'data') else df
            
            f.write(f"Total motif sites: {len(data)}\n")
            f.write(f"Unique consolidated motifs: {len(set(row.get('consolidated_motif_id') for row in data))}\n")
            f.write(f"Lines analyzed: {len(set(row.get('line') for row in data))}\n")
            f.write(f"Genes with motifs: {len(set(row.get('gene_id') for row in data))}\n\n")
            
            # Basic motif frequency
            f.write("## Most Common Motifs\n")
            motif_counts = defaultdict(int)
            motif_consensus_map = {}
            for row in data:
                motif_id = row.get('consolidated_motif_id')
                if motif_id:
                    motif_counts[motif_id] += 1
                    if motif_id not in motif_consensus_map:
                        motif_consensus_map[motif_id] = row.get('motif_consensus', '')
            
            # Sort by count and take top 20
            sorted_motifs = sorted(motif_counts.items(), key=lambda x: x[1], reverse=True)[:20]
            for motif, count in sorted_motifs:
                consensus = motif_consensus_map.get(motif, '')
                f.write(f"{motif}\t{consensus}\t{count} sites\n")

def analyze_motif_cluster(motif_id, sites):
    """
    Analyze a single motif cluster to see if sequences are appropriately grouped
    """
    sequences = [site['sequence'] for site in sites]
    streme_patterns = set(site['original_streme_consensus'] for site in sites)
    
    # Calculate sequence diversity metrics
    unique_sequences = set(sequences)
    seq_lengths = [len(seq) for seq in sequences]
    
    # Check base composition variability
    base_compositions = []
    for seq in sequences[:50]:  # Sample first 50 to avoid too much computation
        total_len = len(seq)
        if total_len > 0:
            comp = {
                'A': seq.count('A') / total_len,
                'T': seq.count('T') / total_len,
                'G': seq.count('G') / total_len,
                'C': seq.count('C') / total_len
            }
            base_compositions.append(comp)
    
    # Calculate composition variance (simple measure of diversity)
    if base_compositions:
        avg_A = sum(comp['A'] for comp in base_compositions) / len(base_compositions)
        avg_T = sum(comp['T'] for comp in base_compositions) / len(base_compositions)
        avg_G = sum(comp['G'] for comp in base_compositions) / len(base_compositions)
        avg_C = sum(comp['C'] for comp in base_compositions) / len(base_compositions)
        
        # Calculate variance for each base
        var_A = sum((comp['A'] - avg_A)**2 for comp in base_compositions) / len(base_compositions)
        var_T = sum((comp['T'] - avg_T)**2 for comp in base_compositions) / len(base_compositions)
        var_G = sum((comp['G'] - avg_G)**2 for comp in base_compositions) / len(base_compositions)
        var_C = sum((comp['C'] - avg_C)**2 for comp in base_compositions) / len(base_compositions)
        
        total_variance = var_A + var_T + var_G + var_C
    else:
        total_variance = 0
    
    return {
        'total_sites': len(sites),
        'unique_sequences': len(unique_sequences),
        'sequence_diversity': len(unique_sequences) / len(sites) if sites else 0,
        'streme_patterns_merged': len(streme_patterns),
        'streme_patterns': list(streme_patterns),
        'length_range': (min(seq_lengths), max(seq_lengths)) if seq_lengths else (0, 0),
        'composition_variance': total_variance,
        'sample_sequences': sequences[:5],  # First 5 sequences as examples
        'consensus': sites[0]['motif_consensus'] if sites else '',
        'cluster_size': sites[0]['cluster_size'] if sites else 0
    }

def run_validation(tsv_file):
    """
    Validate motif consolidation by analyzing the consolidated_streme_sites.tsv output.
    """
    print("=== MOTIF CONSOLIDATION VALIDATION ===")
    print(f"Analyzing: {tsv_file}")
    print()
    
    # Read the consolidated data
    motif_clusters = defaultdict(list)
    
    try:
        with open(tsv_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                motif_id = row['consolidated_motif_id']
                motif_clusters[motif_id].append(row)
    except Exception as e:
        print(f"Error reading file: {e}")
        return False
    
    print(f"Found {len(motif_clusters)} consolidated motifs")
    print()
    
    # Analyze each cluster
    suspicious_clusters = []
    
    for motif_id, sites in motif_clusters.items():
        analysis = analyze_motif_cluster(motif_id, sites)
        
        # Flag potentially problematic clusters
        is_suspicious = False
        issues = []
        
        # High sequence diversity might indicate unrelated sequences
        if analysis['sequence_diversity'] > 0.8 and analysis['total_sites'] > 20:
            is_suspicious = True
            issues.append(f"High sequence diversity ({analysis['sequence_diversity']:.2f})")
        
        # Multiple very different STREME patterns merged
        if analysis['streme_patterns_merged'] > 3:
            is_suspicious = True
            issues.append(f"Many STREME patterns merged ({analysis['streme_patterns_merged']})")
        
        # Very high composition variance
        if analysis['composition_variance'] > 0.3:
            is_suspicious = True
            issues.append(f"High composition variance ({analysis['composition_variance']:.3f})")
        
        # Large length range might indicate different motif types
        length_diff = analysis['length_range'][1] - analysis['length_range'][0]
        if length_diff > 5:
            is_suspicious = True
            issues.append(f"Large length range ({analysis['length_range']})")
        
        if is_suspicious:
            suspicious_clusters.append((motif_id, analysis, issues))
    
    # Report results
    print(f"=== VALIDATION RESULTS ===")
    print(f"Total motifs analyzed: {len(motif_clusters)}")
    print(f"Potentially problematic clusters: {len(suspicious_clusters)}")
    print()
    
    if suspicious_clusters:
        print("=== SUSPICIOUS CLUSTERS (May need review) ===")
        
        # Sort by total sites (most common first)
        suspicious_clusters.sort(key=lambda x: x[1]['total_sites'], reverse=True)
        
        for motif_id, analysis, issues in suspicious_clusters[:10]:  # Show top 10
            print(f"\n🚨 {motif_id} ({analysis['total_sites']} sites)")
            print(f"   True consensus: {analysis['consensus']}")
            print(f"   STREME patterns merged: {analysis['streme_patterns']}")
            print(f"   Issues: {', '.join(issues)}")
            print(f"   Sample sequences:")
            for i, seq in enumerate(analysis['sample_sequences'], 1):
                print(f"     {i}. {seq}")
            
            # Recommendation
            if analysis['composition_variance'] > 0.4:
                print("   💡 Recommendation: Very high variance - likely unrelated sequences")
            elif analysis['streme_patterns_merged'] > 5:
                print("   💡 Recommendation: Too many patterns merged - consider stricter threshold")
            else:
                print("   💡 Recommendation: Review manually - may be legitimate variation")
    
    else:
        print("✅ No obviously suspicious clusters found!")
        print("   The consolidation appears to be working appropriately.")
    
    # Summary statistics
    print(f"\n=== SUMMARY STATISTICS ===")
    
    # Overall cluster size distribution
    cluster_sizes = [len(sites) for sites in motif_clusters.values()]
    single_site_clusters = sum(1 for size in cluster_sizes if size == 1)
    small_clusters = sum(1 for size in cluster_sizes if 2 <= size <= 10)
    medium_clusters = sum(1 for size in cluster_sizes if 11 <= size <= 100)
    large_clusters = sum(1 for size in cluster_sizes if size > 100)
    
    print(f"Single-site motifs: {single_site_clusters}")
    print(f"Small clusters (2-10 sites): {small_clusters}")
    print(f"Medium clusters (11-100 sites): {medium_clusters}")
    print(f"Large clusters (>100 sites): {large_clusters}")
    
    # Merged STREME patterns statistics
    streme_merges = [analysis['streme_patterns_merged'] for _, sites in motif_clusters.items() 
                     for analysis in [analyze_motif_cluster(_, sites)]]
    avg_merges = sum(streme_merges) / len(streme_merges) if streme_merges else 0
    max_merges = max(streme_merges) if streme_merges else 0
    
    print(f"Average STREME patterns per cluster: {avg_merges:.1f}")
    print(f"Maximum STREME patterns merged: {max_merges}")
    
    if len(suspicious_clusters) / len(motif_clusters) < 0.05:
        print(f"\n✅ OVERALL ASSESSMENT: Good consolidation quality")
        print(f"   Less than 5% of clusters flagged as suspicious")
        return True
    elif len(suspicious_clusters) / len(motif_clusters) < 0.15:
        print(f"\n⚠️  OVERALL ASSESSMENT: Moderate consolidation quality")
        print(f"   Some clusters may need review")
        return True
    else:
        print(f"\n🚨 OVERALL ASSESSMENT: Many suspicious clusters")
        print(f"   Consider adjusting consolidation parameters")
        return False

def run_consolidation(args):
    """
    Run the motif consolidation process
    """
    # Validate input directory
    results_dir = Path(args.streme_results_dir)
    if not results_dir.exists():
        print(f"Error: Directory {results_dir} does not exist")
        sys.exit(1)
    
    print("=== STREME Sites Consolidator ===")
    print(f"Input directory: {results_dir}")
    print(f"Similarity threshold: {args.threshold}")
    print(f"Overlap merging: {'Enabled' if args.merge_overlaps else 'Disabled'}")
    if args.merge_overlaps:
        print(f"Overlap threshold: {args.overlap_threshold}")
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
    
    # Merge overlapping motifs (removes STREME sliding window artifacts)
    if args.merge_overlaps:
        consolidated_data = merge_overlapping_motifs(consolidated_data, args.overlap_threshold)
    
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
    if HAS_PANDAS and hasattr(df, 'nunique'):
        print(f"Total consolidated motifs: {df['consolidated_motif_id'].nunique()}")
    else:
        print(f"Total consolidated motifs: {df.nunique('consolidated_motif_id')}")
    print(f"Total sites: {len(df)}")
    
    return output_file

def main():
    parser = argparse.ArgumentParser(
        description='STREME Sites Consolidator - Consolidate and validate motif data',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    
    # Create subparsers for different commands
    subparsers = parser.add_subparsers(dest='command', help='Available commands')
    
    # Consolidate command
    consolidate_parser = subparsers.add_parser(
        'consolidate',
        help='Consolidate STREME sites.tsv files across multiple lines',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Basic consolidation
  %(prog)s consolidate /path/to/streme/results/ --output consolidated_sites
  
  # With custom similarity threshold
  %(prog)s consolidate /path/to/streme/results/ --threshold 0.8 --output results/
  
  # Process specific lines only
  %(prog)s consolidate /path/to/streme/results/ --lines IM502,IM664 --output targeted/

Input Directory Structure:
  streme_results/
  ├── streme_IM502/sites.tsv
  ├── streme_IM664/sites.tsv
  └── streme_IM767/sites.tsv
        """
    )
    
    consolidate_parser.add_argument(
        'streme_results_dir',
        help='Directory containing STREME output folders with sites.tsv files'
    )
    
    consolidate_parser.add_argument(
        '--output', '-o',
        default='consolidated_streme_sites',
        help='Output file prefix (default: consolidated_streme_sites)'
    )
    
    consolidate_parser.add_argument(
        '--threshold', '-t',
        type=float,
        default=0.75,
        help='Similarity threshold for motif consolidation (default: 0.75)'
    )
    
    consolidate_parser.add_argument(
        '--lines', '-l',
        help='Comma-separated list of specific lines to process (e.g., IM502,IM664)'
    )
    
    consolidate_parser.add_argument(
        '--merge-overlaps', '-m',
        action='store_true',
        default=True,
        help='Merge overlapping motif hits within the same gene (default: True)'
    )
    
    consolidate_parser.add_argument(
        '--no-merge-overlaps',
        action='store_false',
        dest='merge_overlaps',
        help='Disable merging of overlapping motif hits'
    )
    
    consolidate_parser.add_argument(
        '--overlap-threshold',
        type=float,
        default=0.5,
        help='Minimum overlap fraction to merge motifs (default: 0.5)'
    )
    
    consolidate_parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output'
    )
    
    # Validate command
    validate_parser = subparsers.add_parser(
        'validate',
        help='Validate consolidated motif results for quality',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Validate consolidation results
  %(prog)s validate consolidated_streme_sites.tsv
  
  # Validate and show detailed analysis
  %(prog)s validate results/consolidated_streme_sites.tsv --verbose
        """
    )
    
    validate_parser.add_argument(
        'tsv_file',
        help='Consolidated TSV file to validate (output from consolidate command)'
    )
    
    validate_parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Show detailed analysis for all clusters'
    )
    
    args = parser.parse_args()
    
    # Handle different commands
    if not args.command:
        parser.print_help()
        sys.exit(1)
    
    if args.command == 'consolidate':
        output_file = run_consolidation(args)
        
        # Ask if user wants to validate
        print(f"\n=== VALIDATION OPTION ===")
        response = input(f"Would you like to validate the results in {output_file}? (y/n): ").strip().lower()
        if response in ['y', 'yes']:
            print()
            run_validation(output_file)
    
    elif args.command == 'validate':
        if not Path(args.tsv_file).exists():
            print(f"Error: File {args.tsv_file} does not exist")
            sys.exit(1)
        
        success = run_validation(args.tsv_file)
        sys.exit(0 if success else 1)

if __name__ == "__main__":
    main()
