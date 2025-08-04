#!/usr/bin/env python3
"""
Script to group sequences based on motif presence for pair-wise analysis.
"""

import sys
from pathlib import Path
from Bio import SeqIO
import re

def read_motifs(motif_file):
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

def group_sequences_by_pairs(sequences, motifs):
    """
    Group sequences based on motif presence for pair-wise analysis.
    """
    groups = {}
    
    # Initialize groups for each motif
    for motif in motifs:
        groups[motif['id']] = []
    
    # Add an 'orphan' group for sequences without clear motif matches
    groups['orphan'] = []
    
    # Process sequences in pairs
    seq_list = list(sequences)
    
    for i in range(0, len(seq_list), 2):
        if i + 1 < len(seq_list):
            # Process pair of sequences
            seq1 = seq_list[i]
            seq2 = seq_list[i + 1]
            
            pair_name = f"{seq1.id}_and_{seq2.id}"
            pair_sequences = [seq1, seq2]
            
            # Check which motifs are represented in this pair
            motif_found = False
            for motif in motifs:
                # Simple heuristic: check if sequences have similar characteristics
                # that might indicate shared motifs (you can refine this logic)
                if has_potential_motif_similarity(seq1, seq2):
                    groups[motif['id']].append({
                        'pair_name': pair_name,
                        'sequences': pair_sequences,
                        'motif': motif['id']
                    })
                    motif_found = True
                    break
            
            if not motif_found:
                groups['orphan'].append({
                    'pair_name': pair_name,
                    'sequences': pair_sequences,
                    'motif': 'none'
                })
        else:
            # Handle odd number of sequences - last sequence becomes orphan
            seq = seq_list[i]
            groups['orphan'].append({
                'pair_name': f"{seq.id}_single",
                'sequences': [seq],
                'motif': 'none'
            })
    
    return groups

def has_potential_motif_similarity(seq1, seq2):
    """
    Simple heuristic to check if two sequences might share motifs.
    You can replace this with more sophisticated motif detection.
    """
    # Check GC content similarity
    gc1 = (seq1.seq.count('G') + seq1.seq.count('C')) / len(seq1.seq)
    gc2 = (seq2.seq.count('G') + seq2.seq.count('C')) / len(seq2.seq)
    
    gc_diff = abs(gc1 - gc2)
    
    # Check length similarity
    len_ratio = min(len(seq1.seq), len(seq2.seq)) / max(len(seq1.seq), len(seq2.seq))
    
    # Simple similarity criteria
    return gc_diff < 0.2 and len_ratio > 0.8

def write_grouped_sequences(groups, output_file):
    """
    Write grouped sequences to output file.
    """
    with open(output_file, 'w') as f:
        f.write("# Grouped Sequences by Motif Presence\n")
        f.write("# Generated from pair-wise sequence analysis\n")
        f.write("# " + "="*60 + "\n\n")
        
        for motif_id, group_data in groups.items():
            if group_data:  # Only write groups that have sequences
                f.write(f"## Motif Group: {motif_id}\n")
                f.write(f"Number of pairs/sequences: {len(group_data)}\n\n")
                
                for item in group_data:
                    f.write(f"Pair: {item['pair_name']}\n")
                    f.write(f"Associated Motif: {item['motif']}\n")
                    
                    for seq in item['sequences']:
                        f.write(f"  Sequence ID: {seq.id}\n")
                        f.write(f"  Length: {len(seq.seq)} bp\n")
                        f.write(f"  Description: {seq.description}\n")
                    
                    f.write("-" * 40 + "\n")
                
                f.write("\n")

def main():
    if len(sys.argv) != 4:
        print("Usage: python group_sequences_by_motifs.py <sequences_file> <motifs_file> <output_file>")
        sys.exit(1)
    
    sequences_file = sys.argv[1]
    motifs_file = sys.argv[2]
    output_file = sys.argv[3]
    
    # Check input files exist
    for file_path in [sequences_file, motifs_file]:
        if not Path(file_path).exists():
            print(f"Error: Input file {file_path} does not exist")
            sys.exit(1)
    
    print(f"Reading sequences from: {sequences_file}")
    print(f"Reading motifs from: {motifs_file}")
    
    # Read sequences
    try:
        sequences = list(SeqIO.parse(sequences_file, "fasta"))
        print(f"Loaded {len(sequences)} sequences")
    except Exception as e:
        print(f"Error reading sequences: {e}")
        sys.exit(1)
    
    # Read motifs
    try:
        motifs = read_motifs(motifs_file)
        print(f"Loaded {len(motifs)} motifs")
    except Exception as e:
        print(f"Error reading motifs: {e}")
        sys.exit(1)
    
    # Group sequences by pairs and motifs
    groups = group_sequences_by_pairs(sequences, motifs)
    
    # Count total assignments
    total_assignments = sum(len(group_data) for group_data in groups.values())
    print(f"Created {len(groups)} groups with {total_assignments} total assignments")
    
    # Write output
    write_grouped_sequences(groups, output_file)
    print(f"Grouped sequences written to: {output_file}")

if __name__ == "__main__":
    main()
