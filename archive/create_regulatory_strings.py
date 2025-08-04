#!/usr/bin/env python3
"""
Advanced regulatory string generator for deep learning applications.
Creates enriched sequences that capture motif presence, order, and spacing.
"""

import sys
import json
from pathlib import Path
import numpy as np

def load_compiled_data(json_file):
    """
    Load the compiled motif data from JSON file.
    """
    with open(json_file, 'r') as f:
        data = json.load(f)
    return data['all_motifs'], data['line_motifs'], data['motif_positions']

def create_enhanced_regulatory_strings(all_motifs, line_motifs, motif_positions, output_file, mode='consensus'):
    """
    Create enhanced regulatory strings with different encoding modes.
    
    Modes:
    - 'consensus': Use consensus sequences
    - 'positional': Include relative positions
    - 'spaced': Include spacing information
    - 'encoded': Use simplified motif codes
    """
    
    # Create motif encoding dictionary
    unique_consensuses = set()
    for motif_data in all_motifs.values():
        if motif_data['consensus']:
            unique_consensuses.add(motif_data['consensus'])
    
    motif_codes = {}
    for i, consensus in enumerate(sorted(unique_consensuses), 1):
        motif_codes[consensus] = f"M{i:03d}"
    
    with open(output_file, 'w') as f:
        f.write("# Enhanced Regulatory Strings for Deep Learning\n")
        f.write(f"# Encoding mode: {mode}\n")
        f.write("# " + "="*60 + "\n\n")
        
        for line_name, motifs in line_motifs.items():
            # Collect motif instances with positions
            motif_instances = []
            
            for i, motif in enumerate(motifs):
                motif_key = f"{line_name}_motif_{i+1}_{motif['motif_id']}"
                
                if motif_key in motif_positions:
                    for pos_info in motif_positions[motif_key]:
                        motif_instances.append({
                            'position': pos_info['position'],
                            'consensus': motif['consensus'],
                            'sequence': pos_info['sequence'],
                            'strand': pos_info['strand'],
                            'evalue': motif['evalue'],
                            'width': motif['width']
                        })
            
            # Sort by position (5' to 3')
            motif_instances.sort(key=lambda x: x['position'])
            
            if not motif_instances:
                continue
            
            # Generate regulatory string based on mode
            if mode == 'consensus':
                reg_string = create_consensus_string(motif_instances)
            elif mode == 'positional':
                reg_string = create_positional_string(motif_instances)
            elif mode == 'spaced':
                reg_string = create_spaced_string(motif_instances)
            elif mode == 'encoded':
                reg_string = create_encoded_string(motif_instances, motif_codes)
            else:
                reg_string = create_consensus_string(motif_instances)
            
            # Write FASTA entry
            header = f">{line_name}_regulatory"
            header += f" |mode={mode}"
            header += f" |motifs={len(motif_instances)}"
            header += f" |length={len(reg_string)}"
            
            # Add metadata
            consensuses = [m['consensus'] for m in motif_instances]
            positions = [str(m['position']) for m in motif_instances]
            
            header += f" |consensuses={','.join(consensuses)}"
            header += f" |positions={','.join(positions)}"
            
            f.write(header + "\n")
            f.write(reg_string + "\n")

def create_consensus_string(motif_instances):
    """
    Create string using consensus sequences separated by spacers.
    """
    sequences = []
    for instance in motif_instances:
        if instance['strand'] == '-':
            # Add indicator for reverse strand
            sequences.append(f"[{instance['consensus']}]")
        else:
            sequences.append(instance['consensus'])
    
    return "NNN".join(sequences)  # Use NNN as spacer

def create_positional_string(motif_instances):
    """
    Create string that includes relative position information.
    """
    if not motif_instances:
        return ""
    
    # Calculate relative positions (normalized to 0-1000 scale)
    min_pos = min(m['position'] for m in motif_instances)
    max_pos = max(m['position'] for m in motif_instances)
    pos_range = max_pos - min_pos if max_pos > min_pos else 1
    
    sequences = []
    for instance in motif_instances:
        rel_pos = int(((instance['position'] - min_pos) / pos_range) * 1000)
        pos_code = f"P{rel_pos:03d}"
        
        if instance['strand'] == '-':
            sequences.append(f"{pos_code}[{instance['consensus']}]")
        else:
            sequences.append(f"{pos_code}{instance['consensus']}")
    
    return "NN".join(sequences)

def create_spaced_string(motif_instances):
    """
    Create string that represents actual spacing between motifs.
    """
    if not motif_instances:
        return ""
    
    if len(motif_instances) == 1:
        return motif_instances[0]['consensus']
    
    sequences = [motif_instances[0]['consensus']]
    
    for i in range(1, len(motif_instances)):
        prev_end = motif_instances[i-1]['position'] + motif_instances[i-1]['width']
        curr_start = motif_instances[i]['position']
        spacing = max(0, curr_start - prev_end)
        
        # Represent spacing with proportional Ns (max 50 for readability)
        spacer_length = min(spacing, 50)
        spacer = "N" * spacer_length
        
        sequences.append(spacer)
        sequences.append(motif_instances[i]['consensus'])
    
    return "".join(sequences)

def create_encoded_string(motif_instances, motif_codes):
    """
    Create string using simplified motif codes.
    """
    codes = []
    for instance in motif_instances:
        consensus = instance['consensus']
        if consensus in motif_codes:
            code = motif_codes[consensus]
            if instance['strand'] == '-':
                code = f"[{code}]"
            codes.append(code)
    
    return "_".join(codes)

def create_motif_vocabulary(all_motifs, output_file):
    """
    Create a vocabulary file of all unique motifs for ML applications.
    """
    unique_motifs = {}
    
    for motif_key, motif_data in all_motifs.items():
        consensus = motif_data['consensus']
        if consensus and consensus not in unique_motifs:
            unique_motifs[consensus] = {
                'width': motif_data['width'],
                'examples': []
            }
        
        if consensus:
            unique_motifs[consensus]['examples'].append({
                'line': motif_data['line'],
                'evalue': motif_data['evalue'],
                'sites': motif_data['sites']
            })
    
    with open(output_file, 'w') as f:
        f.write("# Motif Vocabulary for Deep Learning\n")
        f.write("# Format: Consensus\tWidth\tFrequency\tBest_Evalue\tLines\n")
        f.write("# " + "="*60 + "\n")
        
        for i, (consensus, data) in enumerate(sorted(unique_motifs.items()), 1):
            examples = data['examples']
            frequency = len(examples)
            best_evalue = min(ex['evalue'] for ex in examples)
            lines = set(ex['line'] for ex in examples)
            
            f.write(f"M{i:03d}\t{consensus}\t{data['width']}\t{frequency}\t{best_evalue:.2e}\t{','.join(sorted(lines))}\n")

def create_training_matrix(all_motifs, line_motifs, output_file):
    """
    Create a binary matrix suitable for machine learning training.
    """
    # Get unique consensuses
    unique_consensuses = sorted(set(
        motif_data['consensus'] 
        for motif_data in all_motifs.values() 
        if motif_data['consensus']
    ))
    
    lines = sorted(line_motifs.keys())
    
    with open(output_file, 'w') as f:
        # Header
        f.write("Line\t" + "\t".join(unique_consensuses) + "\n")
        
        for line in lines:
            # Find which motifs are present in this line
            line_consensuses = set()
            
            for motif in line_motifs[line]:
                if motif['consensus']:
                    line_consensuses.add(motif['consensus'])
            
            # Create binary vector
            presence_vector = []
            for consensus in unique_consensuses:
                presence_vector.append("1" if consensus in line_consensuses else "0")
            
            f.write(f"{line}\t" + "\t".join(presence_vector) + "\n")

def main():
    if len(sys.argv) != 3:
        print("Usage: python create_regulatory_strings.py <compiled_data.json> <output_prefix>")
        print("Example: python create_regulatory_strings.py results_raw_data.json ml_ready")
        sys.exit(1)
    
    json_file = sys.argv[1]
    output_prefix = sys.argv[2]
    
    if not Path(json_file).exists():
        print(f"Error: JSON file {json_file} does not exist")
        print("Run compile_all_motifs.py first to generate the compiled data")
        sys.exit(1)
    
    print("Loading compiled motif data...")
    all_motifs, line_motifs, motif_positions = load_compiled_data(json_file)
    
    print(f"Loaded data for {len(line_motifs)} lines with {len(all_motifs)} total motifs")
    
    # Create different encoding modes
    modes = ['consensus', 'positional', 'spaced', 'encoded']
    
    for mode in modes:
        print(f"Creating regulatory strings with {mode} encoding...")
        output_file = f"{output_prefix}_regulatory_{mode}.fasta"
        create_enhanced_regulatory_strings(all_motifs, line_motifs, motif_positions, output_file, mode)
    
    # Create supporting files for ML
    print("Creating motif vocabulary...")
    create_motif_vocabulary(all_motifs, f"{output_prefix}_motif_vocabulary.txt")
    
    print("Creating training matrix...")
    create_training_matrix(all_motifs, line_motifs, f"{output_prefix}_training_matrix.txt")
    
    print(f"Deep learning files created:")
    for mode in modes:
        print(f"  - {output_prefix}_regulatory_{mode}.fasta")
    print(f"  - {output_prefix}_motif_vocabulary.txt")
    print(f"  - {output_prefix}_training_matrix.txt")

if __name__ == "__main__":
    main()
