#!/usr/bin/env python3
"""
Simple script to group sequences based on motif presence for pair-wise analysis.
No external dependencies required.
"""

import sys
from pathlib import Path

def read_fasta_simple(fasta_file):
    """
    Simple FASTA reader without BioPython dependency.
    """
    sequences = []
    current_id = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                # Save previous sequence
                if current_id is not None:
                    sequences.append({
                        'id': current_id,
                        'description': current_id,
                        'sequence': ''.join(current_seq)
                    })
                
                # Start new sequence
                current_id = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)
        
        # Save last sequence
        if current_id is not None:
            sequences.append({
                'id': current_id,
                'description': current_id,
                'sequence': ''.join(current_seq)
            })
    
    return sequences

def read_motifs_simple(motif_file):
    """
    Read consolidated motifs from file.
    """
    motifs = []
    
    with open(motif_file, 'r') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            
            parts = line.strip().split('\t')
            if len(parts) >= 2:
                motifs.append({
                    'id': parts[0],
                    'name': parts[1] if len(parts) > 1 else parts[0],
                    'evalue': parts[2] if len(parts) > 2 else 'N/A',
                    'sites': parts[3] if len(parts) > 3 else 'N/A',
                    'width': parts[4] if len(parts) > 4 else 'N/A'
                })
    
    return motifs

def calculate_gc_content(sequence):
    """
    Calculate GC content of a sequence.
    """
    if not sequence:
        return 0
    
    gc_count = sequence.upper().count('G') + sequence.upper().count('C')
    return gc_count / len(sequence)

def has_sequence_similarity(seq1, seq2):
    """
    Check if two sequences have similar characteristics.
    """
    # GC content similarity
    gc1 = calculate_gc_content(seq1['sequence'])
    gc2 = calculate_gc_content(seq2['sequence'])
    gc_diff = abs(gc1 - gc2)
    
    # Length similarity
    len1, len2 = len(seq1['sequence']), len(seq2['sequence'])
    len_ratio = min(len1, len2) / max(len1, len2) if max(len1, len2) > 0 else 0
    
    # Simple similarity criteria
    return gc_diff < 0.2 and len_ratio > 0.8

def group_sequences_by_pairs_simple(sequences, motifs):
    """
    Group sequences based on pair-wise analysis.
    """
    groups = {}
    
    # Initialize groups
    for motif in motifs:
        groups[motif['id']] = []
    groups['similar_pairs'] = []
    groups['orphans'] = []
    
    # Process sequences in pairs
    for i in range(0, len(sequences), 2):
        if i + 1 < len(sequences):
            # Process pair
            seq1 = sequences[i]
            seq2 = sequences[i + 1]
            
            pair_name = f"pair_{i//2 + 1}_{seq1['id']}_and_{seq2['id']}"
            
            if has_sequence_similarity(seq1, seq2):
                groups['similar_pairs'].append({
                    'pair_name': pair_name,
                    'sequences': [seq1, seq2],
                    'gc1': calculate_gc_content(seq1['sequence']),
                    'gc2': calculate_gc_content(seq2['sequence']),
                    'len1': len(seq1['sequence']),
                    'len2': len(seq2['sequence'])
                })
            else:
                # Add to individual orphans
                groups['orphans'].extend([
                    {
                        'id': seq1['id'],
                        'sequence': seq1,
                        'reason': 'dissimilar_pair'
                    },
                    {
                        'id': seq2['id'],
                        'sequence': seq2,
                        'reason': 'dissimilar_pair'
                    }
                ])
        else:
            # Single sequence left
            seq = sequences[i]
            groups['orphans'].append({
                'id': seq['id'],
                'sequence': seq,
                'reason': 'single_remaining'
            })
    
    return groups

def write_grouped_output_simple(groups, output_file):
    """
    Write grouped sequences to output file.
    """
    with open(output_file, 'w') as f:
        f.write("# Sequence Grouping Results - Pair-wise Analysis\n")
        f.write("# STREME-based promoter sequence consolidation\n")
        f.write("# " + "="*60 + "\n\n")
        
        # Write similar pairs
        if groups['similar_pairs']:
            f.write("## SIMILAR PAIRS\n")
            f.write(f"Found {len(groups['similar_pairs'])} similar sequence pairs\n\n")
            
            for pair in groups['similar_pairs']:
                f.write(f"Pair: {pair['pair_name']}\n")
                f.write(f"Sequence 1: {pair['sequences'][0]['id']}\n")
                f.write(f"  Length: {pair['len1']} bp, GC: {pair['gc1']:.3f}\n")
                f.write(f"Sequence 2: {pair['sequences'][1]['id']}\n")
                f.write(f"  Length: {pair['len2']} bp, GC: {pair['gc2']:.3f}\n")
                f.write("-" * 50 + "\n")
            f.write("\n")
        
        # Write orphans
        if groups['orphans']:
            f.write("## ORPHAN SEQUENCES\n")
            f.write(f"Found {len(groups['orphans'])} orphan sequences\n\n")
            
            for orphan in groups['orphans']:
                seq = orphan['sequence']
                f.write(f"Sequence: {orphan['id']}\n")
                f.write(f"  Length: {len(seq['sequence'])} bp\n")
                f.write(f"  GC content: {calculate_gc_content(seq['sequence']):.3f}\n")
                f.write(f"  Reason: {orphan['reason']}\n")
                f.write("-" * 30 + "\n")
        
        # Write summary statistics
        f.write("\n## SUMMARY\n")
        f.write(f"Total similar pairs: {len(groups['similar_pairs'])}\n")
        f.write(f"Total orphan sequences: {len(groups['orphans'])}\n")
        f.write(f"Total sequences processed: {len(groups['similar_pairs']) * 2 + len(groups['orphans'])}\n")

def main():
    if len(sys.argv) != 4:
        print("Usage: python group_sequences_simple.py <sequences_file> <motifs_file> <output_file>")
        sys.exit(1)
    
    sequences_file = sys.argv[1]
    motifs_file = sys.argv[2]
    output_file = sys.argv[3]
    
    # Check input files
    for file_path in [sequences_file, motifs_file]:
        if not Path(file_path).exists():
            print(f"Error: Input file {file_path} does not exist")
            sys.exit(1)
    
    print(f"Reading sequences from: {sequences_file}")
    print(f"Reading motifs from: {motifs_file}")
    
    # Read data
    sequences = read_fasta_simple(sequences_file)
    motifs = read_motifs_simple(motifs_file)
    
    print(f"Loaded {len(sequences)} sequences")
    print(f"Loaded {len(motifs)} motifs")
    
    # Group sequences
    groups = group_sequences_by_pairs_simple(sequences, motifs)
    
    # Write output
    write_grouped_output_simple(groups, output_file)
    print(f"Grouping results written to: {output_file}")

if __name__ == "__main__":
    main()
