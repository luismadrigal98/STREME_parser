#!/usr/bin/env python3
"""
Quick test script to verify motif parsing works correctly
Usage: python test_motif_parsing.py <path_to_streme_directory>
"""

import sys
import os
import re
from pathlib import Path

def parse_streme_file_test(streme_file):
    """
    Test version of STREME parsing function
    """
    motifs = []
    
    try:
        with open(streme_file, 'r') as f:
            content = f.read()
        
        # Show a sample of the file content
        print("First 1000 characters of file:")
        print(content[:1000])
        print("\n" + "="*60 + "\n")
        
        # Look for all MOTIF lines first
        motif_lines = re.findall(r'MOTIF.*', content)
        print(f"Found {len(motif_lines)} MOTIF lines:")
        for i, line in enumerate(motif_lines[:5]):
            print(f"  {i+1}: {line}")
        print("\n" + "="*60 + "\n")
        
        # Try different parsing approaches
        print("Trying original regex pattern...")
        motif_sections = re.findall(r'MOTIF\s+(\d+)-([ATCGRYSWKMBDHVN]+)\s+(STREME-\d+)(.*?)(?=MOTIF|\Z)', content, re.DOTALL)
        print(f"Original pattern found: {len(motif_sections)} sections")
        
        print("Trying simpler pattern...")
        simple_sections = re.findall(r'MOTIF\s+(\S+)(.*?)(?=MOTIF|\Z)', content, re.DOTALL)
        print(f"Simple pattern found: {len(simple_sections)} sections")
        
        if simple_sections:
            print(f"\nFirst simple section motif ID: '{simple_sections[0][0]}'")
            print(f"First 200 chars of first section data:")
            print(simple_sections[0][1][:200])
        
        motifs = []
        
        print(f"\nFound {len(motif_sections)} motif sections")
        
        for i, (motif_number, consensus_seq, streme_id, motif_data) in enumerate(motif_sections[:5]):  # Show first 5
            motif_info = {
                'id': f"{motif_number}-{consensus_seq}",
                'consensus': consensus_seq,
                'streme_id': streme_id,
                'number': int(motif_number),
                'evalue': '',
                'sites': '',
                'width': len(consensus_seq),
            }
            
            # Extract E-value
            evalue_match = re.search(r'E=\s*([0-9.e-]+)', motif_data)
            if evalue_match:
                motif_info['evalue'] = float(evalue_match.group(1))
            
            # Extract sites count
            sites_match = re.search(r'nsites=\s*(\d+)', motif_data)
            if sites_match:
                motif_info['sites'] = int(sites_match.group(1))
            
            motifs.append(motif_info)
            
            print(f"Motif {i+1}:")
            print(f"  ID: {motif_info['id']}")
            print(f"  Consensus: {motif_info['consensus']}")
            print(f"  E-value: {motif_info['evalue']}")
            print(f"  Sites: {motif_info['sites']}")
            print(f"  Width: {motif_info['width']}")
            print()
    
    except Exception as e:
        print(f"Error parsing {streme_file}: {e}")
        return []
    
    return motifs

def main():
    if len(sys.argv) != 2:
        print("Usage: python test_motif_parsing.py <path_to_streme_txt_file>")
        print("Example: python test_motif_parsing.py /path/to/streme_output/streme.txt")
        sys.exit(1)
    
    streme_file = sys.argv[1]
    
    if not Path(streme_file).exists():
        print(f"Error: File {streme_file} does not exist")
        sys.exit(1)
    
    print(f"Testing motif parsing on: {streme_file}")
    print("="*60)
    
    motifs = parse_streme_file_test(streme_file)
    
    print(f"\nSUMMARY:")
    print(f"Total motifs parsed: {len(motifs)}")
    
    if motifs:
        print(f"First motif consensus: {motifs[0]['consensus']}")
        print(f"Last motif consensus: {motifs[-1]['consensus']}")

if __name__ == "__main__":
    main()
