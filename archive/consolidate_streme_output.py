#!/usr/bin/env python3
"""
Script to consolidate STREME output files by processing pairs of lines
for motif discovery results in Mimulus guttatus promoter analysis.
"""

import sys
import re
from pathlib import Path

def parse_streme_output(input_file):
    """
    Parse STREME output file and extract motif information.
    Process pairs of lines to consolidate motif data.
    """
    motifs = []
    
    with open(input_file, 'r') as f:
        lines = f.readlines()
    
    # Process lines in pairs or according to STREME format
    i = 0
    while i < len(lines):
        line = lines[i].strip()
        
        # Look for motif headers
        if line.startswith('MOTIF'):
            motif_info = {
                'id': '',
                'consensus': '',
                'evalue': '',
                'sites': '',
                'width': ''
            }
            
            # Extract motif ID and basic info
            motif_parts = line.split()
            if len(motif_parts) >= 2:
                motif_info['id'] = motif_parts[1]
            
            # Look ahead for additional information
            j = i + 1
            while j < len(lines) and j < i + 10:  # Look ahead up to 10 lines
                next_line = lines[j].strip()
                
                # Extract E-value
                if 'E-value' in next_line or 'E=' in next_line:
                    evalue_match = re.search(r'E[=-]\s*([0-9.e-]+)', next_line)
                    if evalue_match:
                        motif_info['evalue'] = evalue_match.group(1)
                
                # Extract sites count
                sites_match = re.search(r'(\d+)\s+sites', next_line)
                if sites_match:
                    motif_info['sites'] = sites_match.group(1)
                
                # Extract width
                width_match = re.search(r'width\s*=?\s*(\d+)', next_line)
                if width_match:
                    motif_info['width'] = width_match.group(1)
                
                j += 1
            
            motifs.append(motif_info)
            i = j
        else:
            i += 1
    
    return motifs

def consolidate_motifs(motifs):
    """
    Consolidate similar motifs based on pair-wise comparison.
    """
    consolidated = []
    
    for motif in motifs:
        # Add basic filtering criteria
        try:
            evalue = float(motif['evalue']) if motif['evalue'] else 1.0
            sites = int(motif['sites']) if motif['sites'] else 0
            
            # Filter by significance and site count
            if evalue <= 0.05 and sites >= 5:
                consolidated.append(motif)
        except ValueError:
            # Include motifs with missing/invalid numeric data for manual review
            consolidated.append(motif)
    
    return consolidated

def write_consolidated_output(motifs, output_file):
    """
    Write consolidated motifs to output file.
    """
    with open(output_file, 'w') as f:
        f.write("# Consolidated STREME Motifs\n")
        f.write("# Format: MotifID\tE-value\tSites\tWidth\n")
        f.write("# " + "="*50 + "\n")
        
        for i, motif in enumerate(motifs, 1):
            f.write(f"Motif_{i:02d}\t")
            f.write(f"{motif['id']}\t")
            f.write(f"{motif['evalue']}\t")
            f.write(f"{motif['sites']}\t")
            f.write(f"{motif['width']}\n")

def main():
    if len(sys.argv) != 3:
        print("Usage: python consolidate_streme_output.py <input_streme_file> <output_file>")
        sys.exit(1)
    
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    
    if not Path(input_file).exists():
        print(f"Error: Input file {input_file} does not exist")
        sys.exit(1)
    
    print(f"Processing STREME output file: {input_file}")
    
    # Parse STREME output
    motifs = parse_streme_output(input_file)
    print(f"Found {len(motifs)} motifs")
    
    # Consolidate motifs
    consolidated_motifs = consolidate_motifs(motifs)
    print(f"Consolidated to {len(consolidated_motifs)} significant motifs")
    
    # Write output
    write_consolidated_output(consolidated_motifs, output_file)
    print(f"Consolidated motifs written to: {output_file}")

if __name__ == "__main__":
    main()
