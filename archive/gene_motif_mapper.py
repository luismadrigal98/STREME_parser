#!/usr/bin/env python3
"""
Gene Motif Mapper - CLI tool for mapping motifs to individual genes with position scoring
"""

import os
import re
import sys
import argparse
from pathlib import Path
from collections import defaultdict
from Bio import SeqIO

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
                occurrences = int(parts[5])
                lines = parts[6].split(',') if len(parts) > 6 else []
                
                motifs[motif_id] = {
                    'consensus': consensus,
                    'evalue': evalue,
                    'sites': sites,
                    'width': int(width),
                    'occurrences': occurrences,
                    'lines': lines
                }
    
    return motifs

def expand_iupac(sequence):
    """Expand IUPAC codes to regular expression pattern"""
    iupac_map = {
        'A': 'A', 'T': 'T', 'G': 'G', 'C': 'C',
        'R': '[AG]', 'Y': '[CT]', 'S': '[GC]', 'W': '[AT]',
        'K': '[GT]', 'M': '[AC]', 'B': '[CGT]', 'D': '[AGT]',
        'H': '[ACT]', 'V': '[ACG]', 'N': '[ACGT]'
    }
    
    pattern = ''
    for base in sequence.upper():
        pattern += iupac_map.get(base, base)
    
    return pattern

def find_motif_matches(sequence, motif_consensus, allow_mismatch=1):
    """
    Find all matches of motif in sequence with optional mismatches.
    Returns list of (start_pos, end_pos, match_sequence, score) tuples.
    """
    matches = []
    seq_str = str(sequence).upper()
    motif_pattern = expand_iupac(motif_consensus)
    
    # Exact matches using regex
    try:
        regex_pattern = re.compile(motif_pattern)
        for match in regex_pattern.finditer(seq_str):
            matches.append({
                'start': match.start(),
                'end': match.end(),
                'sequence': match.group(),
                'score': 1.0,  # Perfect match
                'type': 'exact'
            })
    except re.error:
        # Fallback if regex is invalid
        pass
    
    # Add fuzzy matching if requested
    if allow_mismatch > 0:
        motif_len = len(motif_consensus)
        for i in range(len(seq_str) - motif_len + 1):
            subseq = seq_str[i:i + motif_len]
            
            # Skip if we already have an exact match here
            if any(match['start'] == i and match['type'] == 'exact' for match in matches):
                continue
            
            # Calculate similarity
            mismatches = 0
            for j, (m_base, s_base) in enumerate(zip(motif_consensus.upper(), subseq)):
                if m_base in 'RYSWKMBDHVN':
                    # IUPAC code - check if it matches
                    possible_bases = expand_iupac(m_base).strip('[]')
                    if s_base not in possible_bases:
                        mismatches += 1
                elif m_base != s_base:
                    mismatches += 1
            
            # Accept if mismatches are within threshold
            if mismatches <= allow_mismatch:
                score = 1.0 - (mismatches / motif_len)
                matches.append({
                    'start': i,
                    'end': i + motif_len,
                    'sequence': subseq,
                    'score': score,
                    'type': 'fuzzy'
                })
    
    return matches

def calculate_position_weight(position, gene_length, region_type='promoter'):
    """
    Calculate position weight based on distance from gene start.
    Higher weights for regulatory important regions.
    """
    if region_type == 'promoter':
        # Promoter region - higher weight closer to TSS
        if position < 500:  # Core promoter
            return 1.0
        elif position < 1000:  # Proximal promoter
            return 0.8
        elif position < 2000:  # Distal promoter
            return 0.6
        else:
            return 0.3
    
    elif region_type == 'gene_body':
        # Gene body - consider introns vs exons
        relative_pos = position / gene_length
        if relative_pos < 0.1:  # 5' UTR region
            return 0.9
        elif relative_pos > 0.9:  # 3' UTR region
            return 0.7
        else:  # Gene body
            return 0.5
    
    return 0.1

def map_motifs_to_genes(fasta_file, motifs, line_name, allow_mismatch=1, min_score=0.7):
    """
    Map all motifs to all genes in a FASTA file.
    Returns detailed gene-motif mapping with positions and scores.
    """
    gene_motif_map = {}
    
    try:
        sequences = SeqIO.parse(fasta_file, 'fasta')
        
        for seq_record in sequences:
            gene_id = seq_record.id
            gene_seq = seq_record.seq
            gene_length = len(gene_seq)
            
            gene_motif_map[gene_id] = {
                'length': gene_length,
                'motifs': [],
                'total_score': 0.0,
                'motif_count': 0
            }
            
            # Check each motif against this gene
            for motif_id, motif_data in motifs.items():
                # Only check motifs that are present in this line
                if line_name not in motif_data['lines']:
                    continue
                
                consensus = motif_data['consensus']
                matches = find_motif_matches(gene_seq, consensus, allow_mismatch)
                
                for match in matches:
                    if match['score'] >= min_score:
                        # Calculate position weight
                        pos_weight = calculate_position_weight(match['start'], gene_length)
                        final_score = match['score'] * pos_weight
                        
                        motif_hit = {
                            'motif_id': motif_id,
                            'consensus': consensus,
                            'start': match['start'],
                            'end': match['end'],
                            'sequence': match['sequence'],
                            'match_score': match['score'],
                            'position_weight': pos_weight,
                            'final_score': final_score,
                            'match_type': match['type']
                        }
                        
                        gene_motif_map[gene_id]['motifs'].append(motif_hit)
                        gene_motif_map[gene_id]['total_score'] += final_score
                        gene_motif_map[gene_id]['motif_count'] += 1
            
            # Sort motifs by position
            gene_motif_map[gene_id]['motifs'].sort(key=lambda x: x['start'])
    
    except Exception as e:
        print(f"Error processing {fasta_file}: {e}")
    
    return gene_motif_map

def write_gene_motif_matrix(all_gene_maps, motifs, output_file):
    """
    Write gene-motif presence/absence matrix for ML/statistical analysis.
    """
    # Collect all genes and motifs
    all_genes = set()
    all_motif_ids = set()
    
    for line_name, gene_map in all_gene_maps.items():
        all_genes.update(gene_map.keys())
        for gene_data in gene_map.values():
            for motif_hit in gene_data['motifs']:
                all_motif_ids.add(motif_hit['motif_id'])
    
    all_genes = sorted(all_genes)
    all_motif_ids = sorted(all_motif_ids)
    
    with open(output_file, 'w') as f:
        # Header
        f.write("Gene\tLine\t" + "\t".join(all_motif_ids) + "\n")
        
        # Data rows
        for gene_id in all_genes:
            for line_name, gene_map in sorted(all_gene_maps.items()):
                if gene_id in gene_map:
                    gene_data = gene_map[gene_id]
                    
                    # Create presence vector
                    presence_vector = []
                    gene_motifs = {hit['motif_id']: hit['final_score'] for hit in gene_data['motifs']}
                    
                    for motif_id in all_motif_ids:
                        score = gene_motifs.get(motif_id, 0.0)
                        presence_vector.append(f"{score:.3f}")
                    
                    f.write(f"{gene_id}\t{line_name}\t" + "\t".join(presence_vector) + "\n")

def write_detailed_gene_motif_report(all_gene_maps, output_file):
    """
    Write detailed report with all motif positions and scores for each gene.
    """
    with open(output_file, 'w') as f:
        f.write("# Detailed Gene-Motif Mapping Report\n")
        f.write("# Gene\tLine\tMotifID\tConsensus\tStart\tEnd\tSequence\tMatchScore\tPositionWeight\tFinalScore\tType\n")
        
        for line_name, gene_map in sorted(all_gene_maps.items()):
            for gene_id, gene_data in sorted(gene_map.items()):
                if gene_data['motifs']:
                    for motif_hit in gene_data['motifs']:
                        f.write(f"{gene_id}\t{line_name}\t{motif_hit['motif_id']}\t"
                               f"{motif_hit['consensus']}\t{motif_hit['start']}\t{motif_hit['end']}\t"
                               f"{motif_hit['sequence']}\t{motif_hit['match_score']:.3f}\t"
                               f"{motif_hit['position_weight']:.3f}\t{motif_hit['final_score']:.3f}\t"
                               f"{motif_hit['match_type']}\n")

def main():
    parser = argparse.ArgumentParser(
        description="Map motifs to individual genes with position-aware scoring",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Map motifs to genes in one line
  %(prog)s motif_catalog.txt /path/to/sequences/ --line IM502
  
  # Process all lines with custom parameters
  %(prog)s motif_catalog.txt /path/to/sequences/ --all-lines --mismatch 2 --min-score 0.6
        """
    )
    
    parser.add_argument(
        'motif_catalog',
        help='Motif catalog file from motif_catalog_builder.py'
    )
    
    parser.add_argument(
        'sequence_dir',
        help='Directory containing FASTA files for each line'
    )
    
    parser.add_argument(
        '--line', '-l',
        help='Specific line to process (e.g., "IM502")'
    )
    
    parser.add_argument(
        '--all-lines', '-a',
        action='store_true',
        help='Process all lines found in sequence directory'
    )
    
    parser.add_argument(
        '--output', '-o',
        default='gene_motif_mapping',
        help='Output file prefix (default: gene_motif_mapping)'
    )
    
    parser.add_argument(
        '--mismatch', '-m',
        type=int,
        default=1,
        help='Allow mismatches in motif matching (default: 1)'
    )
    
    parser.add_argument(
        '--min-score', '-s',
        type=float,
        default=0.7,
        help='Minimum match score to report (default: 0.7)'
    )
    
    parser.add_argument(
        '--verbose', '-v',
        action='store_true',
        help='Enable verbose output'
    )
    
    args = parser.parse_args()
    
    if not args.line and not args.all_lines:
        print("Error: Must specify either --line or --all-lines")
        sys.exit(1)
    
    if not os.path.exists(args.motif_catalog):
        print(f"Error: Motif catalog file {args.motif_catalog} does not exist")
        sys.exit(1)
    
    if not os.path.exists(args.sequence_dir):
        print(f"Error: Sequence directory {args.sequence_dir} does not exist")
        sys.exit(1)
    
    # Read motif catalog
    print(f"Reading motif catalog from: {args.motif_catalog}")
    motifs = read_motif_catalog(args.motif_catalog)
    print(f"Loaded {len(motifs)} motifs")
    
    # Find FASTA files to process
    fasta_files = []
    if args.all_lines:
        for file in os.listdir(args.sequence_dir):
            if file.endswith('.fasta') or file.endswith('.fa'):
                fasta_files.append((file, os.path.join(args.sequence_dir, file)))
    else:
        # Look for specific line
        for file in os.listdir(args.sequence_dir):
            if args.line in file and (file.endswith('.fasta') or file.endswith('.fa')):
                fasta_files.append((file, os.path.join(args.sequence_dir, file)))
        
        if not fasta_files:
            print(f"Error: No FASTA files found for line {args.line}")
            sys.exit(1)
    
    if not fasta_files:
        print("Error: No FASTA files found to process")
        sys.exit(1)
    
    print(f"Processing {len(fasta_files)} FASTA files")
    
    # Process each FASTA file
    all_gene_maps = {}
    
    for file_name, file_path in fasta_files:
        # Extract line name from filename
        line_name = file_name.replace('.fasta', '').replace('.fa', '')
        
        if args.verbose:
            print(f"Processing {line_name}...")
        
        gene_map = map_motifs_to_genes(
            file_path, motifs, line_name, 
            args.mismatch, args.min_score
        )
        
        all_gene_maps[line_name] = gene_map
        
        gene_count = len(gene_map)
        motif_hits = sum(len(data['motifs']) for data in gene_map.values())
        
        if args.verbose:
            print(f"  {gene_count} genes, {motif_hits} motif hits")
    
    # Write output files
    matrix_file = f"{args.output}_matrix.txt"
    detail_file = f"{args.output}_detailed.txt"
    
    print(f"\nWriting output files:")
    print(f"- {matrix_file}")
    print(f"- {detail_file}")
    
    write_gene_motif_matrix(all_gene_maps, motifs, matrix_file)
    write_detailed_gene_motif_report(all_gene_maps, detail_file)
    
    # Summary statistics
    total_genes = sum(len(gene_map) for gene_map in all_gene_maps.values())
    total_hits = sum(
        sum(len(gene_data['motifs']) for gene_data in gene_map.values())
        for gene_map in all_gene_maps.values()
    )
    
    print(f"\n=== SUMMARY ===")
    print(f"Total genes processed: {total_genes}")
    print(f"Total motif hits: {total_hits}")
    print(f"Average hits per gene: {total_hits/total_genes:.2f}")

if __name__ == "__main__":
    main()
