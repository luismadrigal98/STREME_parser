#!/usr/bin/env python3
"""
Compile all motifs from individual STREME runs into comprehensive analysis files.
Creates:
1. Master motif catalog
2. Motif presence matrix per line
3. Positional motif maps
4. Regulatory strings for deep learning
"""

import sys
import os
import re
from pathlib import Path
import json

def parse_streme_output(streme_file):
    """
    Parse STREME output file to extract motifs with positions and sequences.
    """
    motifs = []
    
    with open(streme_file, 'r') as f:
        content = f.read()
    
    # Find motif blocks in STREME output
    motif_blocks = re.findall(r'MOTIF\s+(\S+).*?(?=MOTIF|\Z)', content, re.DOTALL)
    
    for i, block in enumerate(motif_blocks):
        lines = block.strip().split('\n')
        if not lines:
            continue
            
        # Extract motif header info
        header_line = lines[0]
        motif_parts = header_line.split()
        
        motif_info = {
            'motif_id': motif_parts[0] if motif_parts else f'motif_{i+1}',
            'width': 0,
            'sites': 0,
            'evalue': 1.0,
            'consensus': '',
            'pwm': [],
            'positions': []
        }
        
        # Parse motif details from subsequent lines
        for line in lines[1:]:
            line = line.strip()
            
            # Extract E-value
            evalue_match = re.search(r'E-value\s*[:=]\s*([0-9.e-]+)', line)
            if evalue_match:
                motif_info['evalue'] = float(evalue_match.group(1))
            
            # Extract width
            width_match = re.search(r'width\s*[:=]\s*(\d+)', line)
            if width_match:
                motif_info['width'] = int(width_match.group(1))
            
            # Extract sites
            sites_match = re.search(r'(\d+)\s+sites', line)
            if sites_match:
                motif_info['sites'] = int(sites_match.group(1))
            
            # Extract consensus sequence
            if len(line) > 0 and all(c in 'ACGT' for c in line.upper()):
                if len(line) > len(motif_info['consensus']):
                    motif_info['consensus'] = line.upper()
        
        motifs.append(motif_info)
    
    return motifs

def find_motif_positions_in_sequences(motif_consensus, fasta_file, max_mismatches=1):
    """
    Find positions of motifs in the original sequences.
    """
    positions = []
    
    # Read sequences
    sequences = read_fasta_simple(fasta_file)
    
    for seq_record in sequences:
        seq = seq_record['sequence'].upper()
        motif = motif_consensus.upper()
        
        # Find exact matches first
        for i in range(len(seq) - len(motif) + 1):
            subseq = seq[i:i+len(motif)]
            
            # Count mismatches
            mismatches = sum(1 for a, b in zip(subseq, motif) if a != b)
            
            if mismatches <= max_mismatches:
                positions.append({
                    'sequence_id': seq_record['id'],
                    'position': i,
                    'sequence': subseq,
                    'strand': '+',
                    'mismatches': mismatches
                })
        
        # Check reverse complement
        rev_comp = reverse_complement(seq)
        for i in range(len(rev_comp) - len(motif) + 1):
            subseq = rev_comp[i:i+len(motif)]
            
            mismatches = sum(1 for a, b in zip(subseq, motif) if a != b)
            
            if mismatches <= max_mismatches:
                # Convert position back to forward strand
                forward_pos = len(seq) - i - len(motif)
                positions.append({
                    'sequence_id': seq_record['id'],
                    'position': forward_pos,
                    'sequence': subseq,
                    'strand': '-',
                    'mismatches': mismatches
                })
    
    return positions

def read_fasta_simple(fasta_file):
    """
    Simple FASTA reader.
    """
    sequences = []
    current_id = None
    current_seq = []
    
    with open(fasta_file, 'r') as f:
        for line in f:
            line = line.strip()
            if line.startswith('>'):
                if current_id is not None:
                    sequences.append({
                        'id': current_id,
                        'sequence': ''.join(current_seq)
                    })
                current_id = line[1:]
                current_seq = []
            else:
                current_seq.append(line)
        
        if current_id is not None:
            sequences.append({
                'id': current_id,
                'sequence': ''.join(current_seq)
            })
    
    return sequences

def reverse_complement(seq):
    """
    Get reverse complement of DNA sequence.
    """
    complement = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N'}
    return ''.join(complement.get(base, base) for base in reversed(seq))

def compile_all_streme_results(streme_dirs, output_prefix):
    """
    Compile all STREME results from multiple directories.
    """
    all_motifs = {}
    line_motifs = {}
    motif_positions = {}
    
    print(f"Processing {len(streme_dirs)} STREME output directories...")
    
    for streme_dir in streme_dirs:
        streme_path = Path(streme_dir)
        streme_file = streme_path / "streme.txt"
        
        if not streme_file.exists():
            print(f"Warning: {streme_file} not found, skipping...")
            continue
        
        # Extract line name from directory
        line_name = streme_path.name.replace('streme_', '')
        print(f"Processing line: {line_name}")
        
        # Parse motifs from this line
        motifs = parse_streme_output(streme_file)
        line_motifs[line_name] = motifs
        
        # Find corresponding FASTA file
        fasta_file = f"{line_name}.fasta.masked"
        if not Path(fasta_file).exists():
            fasta_file = f"{line_name}.fasta"
        
        # Add motifs to master catalog
        for i, motif in enumerate(motifs):
            motif_key = f"{line_name}_motif_{i+1}_{motif['motif_id']}"
            all_motifs[motif_key] = {
                'line': line_name,
                'local_id': motif['motif_id'],
                'consensus': motif['consensus'],
                'width': motif['width'],
                'sites': motif['sites'],
                'evalue': motif['evalue']
            }
            
            # Find positions if FASTA file exists
            if Path(fasta_file).exists():
                positions = find_motif_positions_in_sequences(motif['consensus'], fasta_file)
                motif_positions[motif_key] = positions
    
    return all_motifs, line_motifs, motif_positions

def create_motif_presence_matrix(all_motifs, line_motifs, output_file):
    """
    Create a matrix showing which motifs are present in which lines.
    """
    lines = list(line_motifs.keys())
    unique_consensuses = {}
    
    # Group motifs by consensus sequence
    for motif_key, motif_data in all_motifs.items():
        consensus = motif_data['consensus']
        if consensus not in unique_consensuses:
            unique_consensuses[consensus] = []
        unique_consensuses[consensus].append(motif_key)
    
    with open(output_file, 'w') as f:
        # Header
        f.write("Consensus_Sequence\tWidth\t" + "\t".join(lines) + "\tTotal_Lines\tMotif_IDs\n")
        
        for consensus, motif_keys in unique_consensuses.items():
            if not consensus:  # Skip empty consensuses
                continue
                
            # Count presence in each line
            line_presence = {}
            width = 0
            
            for line in lines:
                line_presence[line] = 0
                
            for motif_key in motif_keys:
                motif_data = all_motifs[motif_key]
                line = motif_data['line']
                line_presence[line] = 1
                width = motif_data['width']
            
            total_lines = sum(line_presence.values())
            
            # Write row
            f.write(f"{consensus}\t{width}\t")
            f.write("\t".join(str(line_presence[line]) for line in lines))
            f.write(f"\t{total_lines}\t")
            f.write(",".join(motif_keys))
            f.write("\n")

def create_positional_maps(motif_positions, line_motifs, output_file):
    """
    Create positional maps showing motif order for each line.
    """
    with open(output_file, 'w') as f:
        f.write("# Positional Motif Maps (5' to 3' order)\n")
        f.write("# Format: Line -> Position -> Motif\n")
        f.write("# " + "="*60 + "\n\n")
        
        for line_name, motifs in line_motifs.items():
            f.write(f"## Line: {line_name}\n")
            
            # Collect all positions for this line
            line_positions = []
            
            for i, motif in enumerate(motifs):
                motif_key = f"{line_name}_motif_{i+1}_{motif['motif_id']}"
                
                if motif_key in motif_positions:
                    for pos_info in motif_positions[motif_key]:
                        line_positions.append({
                            'position': pos_info['position'],
                            'motif_id': motif['motif_id'],
                            'consensus': motif['consensus'],
                            'strand': pos_info['strand'],
                            'sequence': pos_info['sequence']
                        })
            
            # Sort by position (5' to 3')
            line_positions.sort(key=lambda x: x['position'])
            
            f.write(f"Total motifs found: {len(line_positions)}\n")
            f.write("Position\tMotif_ID\tConsensus\tStrand\tSequence\n")
            
            for pos_data in line_positions:
                f.write(f"{pos_data['position']}\t")
                f.write(f"{pos_data['motif_id']}\t")
                f.write(f"{pos_data['consensus']}\t")
                f.write(f"{pos_data['strand']}\t")
                f.write(f"{pos_data['sequence']}\n")
            
            f.write("\n")

def create_regulatory_strings(motif_positions, line_motifs, output_file):
    """
    Create regulatory strings for deep learning applications.
    """
    with open(output_file, 'w') as f:
        f.write("# Regulatory Strings for Deep Learning\n")
        f.write("# Format: FASTA-like with regulatory motif sequences\n")
        f.write("# " + "="*60 + "\n\n")
        
        for line_name, motifs in line_motifs.items():
            # Collect and sort positions
            line_positions = []
            
            for i, motif in enumerate(motifs):
                motif_key = f"{line_name}_motif_{i+1}_{motif['motif_id']}"
                
                if motif_key in motif_positions:
                    for pos_info in motif_positions[motif_key]:
                        line_positions.append({
                            'position': pos_info['position'],
                            'consensus': motif['consensus'],
                            'sequence': pos_info['sequence']
                        })
            
            # Sort by position
            line_positions.sort(key=lambda x: x['position'])
            
            # Create regulatory string
            regulatory_sequence = ""
            motif_labels = []
            
            for i, pos_data in enumerate(line_positions):
                if i > 0:
                    regulatory_sequence += "N"  # Separator between motifs
                regulatory_sequence += pos_data['consensus']
                motif_labels.append(f"pos{pos_data['position']}")
            
            # Write FASTA entry
            f.write(f">{line_name}_regulatory |motifs={len(line_positions)} |labels={','.join(motif_labels)}\n")
            f.write(f"{regulatory_sequence}\n")

def main():
    if len(sys.argv) < 3:
        print("Usage: python compile_all_motifs.py <output_prefix> <streme_dir1> [streme_dir2] ...")
        print("Example: python compile_all_motifs.py results streme_line1 streme_line2 streme_line3")
        sys.exit(1)
    
    output_prefix = sys.argv[1]
    streme_dirs = sys.argv[2:]
    
    # Check that directories exist
    for streme_dir in streme_dirs:
        if not Path(streme_dir).exists():
            print(f"Error: Directory {streme_dir} does not exist")
            sys.exit(1)
    
    print("Starting motif compilation...")
    
    # Compile all results
    all_motifs, line_motifs, motif_positions = compile_all_streme_results(streme_dirs, output_prefix)
    
    print(f"Found {len(all_motifs)} total motifs across {len(line_motifs)} lines")
    
    # Create output files
    print("Creating motif presence matrix...")
    create_motif_presence_matrix(all_motifs, line_motifs, f"{output_prefix}_motif_presence_matrix.txt")
    
    print("Creating positional maps...")
    create_positional_maps(motif_positions, line_motifs, f"{output_prefix}_positional_maps.txt")
    
    print("Creating regulatory strings...")
    create_regulatory_strings(motif_positions, line_motifs, f"{output_prefix}_regulatory_strings.fasta")
    
    # Save raw data as JSON for further analysis
    print("Saving raw data...")
    with open(f"{output_prefix}_raw_data.json", 'w') as f:
        json.dump({
            'all_motifs': all_motifs,
            'line_motifs': line_motifs,
            'motif_positions': motif_positions
        }, f, indent=2)
    
    print(f"Compilation complete! Output files:")
    print(f"  - {output_prefix}_motif_presence_matrix.txt")
    print(f"  - {output_prefix}_positional_maps.txt")
    print(f"  - {output_prefix}_regulatory_strings.fasta")
    print(f"  - {output_prefix}_raw_data.json")

if __name__ == "__main__":
    main()
