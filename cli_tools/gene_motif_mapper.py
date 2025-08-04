#!/usr/bin/env python3
"""
Gene Motif Mapper CLI - Map motifs to individual genes and score presence/position.

This tool takes a motif catalog and gene sequences to create gene-level
regulatory maps suitable for expression analysis.
"""

import os
import re
import sys
import argparse
from pathlib import Path
from collections import defaultdict

def read_motif_catalog(catalog_file):
    """Read consolidated motif catalog"""
    motifs = {}
    
    with open(catalog_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.strip().split('\t')
            if len(parts) >= 6:
                motif_id = parts[0]
                consensus = parts[1]
                evalue = parts[2]
                sites = parts[3]
                width = parts[4]
                occurrences = parts[5]
                lines = parts[6] if len(parts) > 6 else ''
                
                motifs[motif_id] = {
                    'consensus': consensus,
                    'evalue': evalue,
                    'sites': sites,
                    'width': int(width),
                    'occurrences': int(occurrences),
                    'lines': lines.split(',') if lines else []
                }
    
    return motifs

def read_fasta_sequences(fasta_file):
    """Read FASTA sequences without BioPython dependency"""
    sequences = {}
    current_id = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id:
                    sequences[current_id] = ''.join(current_seq)
                current_id = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)
        
        # Don't forget the last sequence
        if current_id:
            sequences[current_id] = ''.join(current_seq)
    
    return sequences

def expand_iupac_pattern(motif_consensus):
    """
    Convert IUPAC consensus to regex pattern for searching
    """
    iupac_map = {
        'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C',
        'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
        'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
        'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
    }
    
    pattern = ''
    for base in motif_consensus.upper():
        pattern += iupac_map.get(base, base)
    
    return pattern

def find_motif_matches(sequence, motif_consensus, min_score=0.8):
    """
    Find all matches of a motif in a sequence
    Returns list of (position, match_sequence, score) tuples
    """
    matches = []
    pattern = expand_iupac_pattern(motif_consensus)
    
    try:
        # Find exact matches first
        for match in re.finditer(pattern, sequence.upper()):
            matches.append({
                'position': match.start(),
                'end_position': match.end(),
                'matched_sequence': match.group(),
                'score': 1.0,  # Exact match
                'strand': '+'
            })
        
        # Also search reverse complement
        rev_comp = reverse_complement(sequence)
        for match in re.finditer(pattern, rev_comp.upper()):
            # Convert position back to forward strand coordinates
            rev_pos = match.start()
            forward_pos = len(sequence) - match.end()
            matches.append({
                'position': forward_pos,
                'end_position': forward_pos + len(match.group()),
                'matched_sequence': match.group(),
                'score': 1.0,  # Exact match
                'strand': '-'
            })
            
    except re.error:
        # If regex pattern is invalid, skip this motif
        pass
    
    return matches

def reverse_complement(sequence):
    """Generate reverse complement of DNA sequence"""
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    rev_comp = ''
    for base in reversed(sequence.upper()):
        rev_comp += complement.get(base, base)
    return rev_comp

def calculate_positional_scores(matches, sequence_length, tss_position=None):
    """
    Calculate position-based scores for motif matches
    Closer to TSS (if provided) or 5' end gets higher scores
    """
    if not matches:
        return []
    
    scored_matches = []
    
    for match in matches:
        # Base score from match quality
        base_score = match['score']
        
        # Position score - closer to 5' end (or TSS) is better
        if tss_position is not None:
            distance_from_tss = abs(match['position'] - tss_position)
            # Score decreases with distance (max 1000bp window)
            position_score = max(0, 1 - distance_from_tss / 1000)
        else:
            # No TSS info - closer to 5' end is better
            position_score = 1 - (match['position'] / sequence_length)
        
        # Combined score
        final_score = base_score * 0.7 + position_score * 0.3
        
        scored_match = match.copy()
        scored_match['position_score'] = position_score
        scored_match['final_score'] = final_score
        
        scored_matches.append(scored_match)
    
    return scored_matches

def map_motifs_to_genes(gene_sequences, motif_catalog, line_name):
    """
    Map all motifs to all genes for a given line
    """
    gene_motif_map = {}
    
    print(f"Mapping motifs to {len(gene_sequences)} genes...")
    
    for gene_idx, (gene_id, sequence) in enumerate(gene_sequences.items()):
        if gene_idx % 1000 == 0:
            print(f"  Processed {gene_idx} genes...")
        
        gene_motif_map[gene_id] = {
            'sequence_length': len(sequence),
            'motifs': {}
        }
        
        # Check each motif against this gene
        for motif_id, motif_info in motif_catalog.items():
            # Skip motifs not found in this line
            if line_name not in motif_info['lines']:
                continue
            
            matches = find_motif_matches(sequence, motif_info['consensus'])
            scored_matches = calculate_positional_scores(matches, len(sequence))
            
            if scored_matches:
                gene_motif_map[gene_id]['motifs'][motif_id] = {
                    'consensus': motif_info['consensus'],
                    'match_count': len(scored_matches),
                    'best_score': max(m['final_score'] for m in scored_matches),
                    'average_score': sum(m['final_score'] for m in scored_matches) / len(scored_matches),
                    'matches': scored_matches
                }
    
    print(f"  Completed mapping for all {len(gene_sequences)} genes")
    return gene_motif_map

def write_gene_motif_matrix(gene_motif_map, motif_catalog, output_file):
    """
    Write gene x motif matrix suitable for expression analysis
    """
    # Get all genes and motifs
    all_genes = sorted(gene_motif_map.keys())
    all_motifs = sorted(motif_catalog.keys())
    
    with open(output_file, 'w') as f:
        # Header
        f.write("GeneID\t" + "\t".join(all_motifs) + "\n")
        
        # Data rows
        for gene_id in all_genes:
            gene_data = gene_motif_map[gene_id]
            scores = []
            
            for motif_id in all_motifs:
                if motif_id in gene_data['motifs']:
                    # Use best score for this gene-motif combination
                    score = gene_data['motifs'][motif_id]['best_score']
                else:
                    score = 0.0
                scores.append(f"{score:.3f}")
            
            f.write(f"{gene_id}\t" + "\t".join(scores) + "\n")

def write_detailed_gene_motif_report(gene_motif_map, output_file):
    """
    Write detailed report with all motif matches per gene
    """
    with open(output_file, 'w') as f:
        f.write("# Detailed Gene-Motif Mapping Report\n")
        f.write("# GeneID\tMotifID\tConsensus\tMatchCount\tBestScore\tAvgScore\tPositions\n")
        
        for gene_id, gene_data in sorted(gene_motif_map.items()):
            for motif_id, motif_data in gene_data['motifs'].items():
                positions = [str(m['position']) for m in motif_data['matches']]
                f.write(f"{gene_id}\t{motif_id}\t{motif_data['consensus']}\t"
                       f"{motif_data['match_count']}\t{motif_data['best_score']:.3f}\t"
                       f"{motif_data['average_score']:.3f}\t{','.join(positions)}\n")

def main():
    parser = argparse.ArgumentParser(
        description='Map motifs to individual genes for expression analysis',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s motif_catalog.txt sequences.fasta IM502 --output outputs/
  %(prog)s motif_catalog.txt sequences.fasta IM502 --output outputs/ --detailed
        """
    )
    
    parser.add_argument('motif_catalog', 
                       help='Consolidated motif catalog file')
    parser.add_argument('gene_sequences',
                       help='FASTA file with gene sequences')
    parser.add_argument('line_name',
                       help='Line name (e.g., IM502, IM664)')
    parser.add_argument('--output', '-o', default='outputs/',
                       help='Output directory (default: outputs/)')
    parser.add_argument('--detailed', action='store_true',
                       help='Generate detailed report with all match positions')
    parser.add_argument('--min-score', type=float, default=0.8,
                       help='Minimum match score threshold (default: 0.8)')
    
    args = parser.parse_args()
    
    # Validate inputs
    for file_path in [args.motif_catalog, args.gene_sequences]:
        if not os.path.exists(file_path):
            print(f"Error: File {file_path} does not exist")
            sys.exit(1)
    
    # Create output directory
    os.makedirs(args.output, exist_ok=True)
    
    print(f"Gene Motif Mapper")
    print(f"Motif catalog: {args.motif_catalog}")
    print(f"Gene sequences: {args.gene_sequences}")
    print(f"Line name: {args.line_name}")
    print(f"Output directory: {args.output}")
    print("="*60)
    
    # Read inputs
    print("Reading motif catalog...")
    motif_catalog = read_motif_catalog(args.motif_catalog)
    print(f"Loaded {len(motif_catalog)} motifs")
    
    print("Reading gene sequences...")
    gene_sequences = read_fasta_sequences(args.gene_sequences)
    print(f"Loaded {len(gene_sequences)} gene sequences")
    
    # Filter motifs for this line
    line_motifs = {mid: minfo for mid, minfo in motif_catalog.items() 
                   if args.line_name in minfo['lines']}
    print(f"Found {len(line_motifs)} motifs for line {args.line_name}")
    
    # Map motifs to genes
    gene_motif_map = map_motifs_to_genes(gene_sequences, line_motifs, args.line_name)
    
    # Count genes with motifs
    genes_with_motifs = sum(1 for g in gene_motif_map.values() if g['motifs'])
    print(f"Found motifs in {genes_with_motifs}/{len(gene_sequences)} genes")
    
    # Write outputs
    matrix_file = os.path.join(args.output, f"{args.line_name}_gene_motif_matrix.txt")
    write_gene_motif_matrix(gene_motif_map, line_motifs, matrix_file)
    print(f"Gene-motif matrix written to: {matrix_file}")
    
    if args.detailed:
        detail_file = os.path.join(args.output, f"{args.line_name}_detailed_motif_report.txt")
        write_detailed_gene_motif_report(gene_motif_map, detail_file)
        print(f"Detailed report written to: {detail_file}")
    
    print("\nDone! Files are ready for expression analysis.")

if __name__ == "__main__":
    main()
