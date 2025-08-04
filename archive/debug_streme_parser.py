#!/usr/bin/env python3
"""
Debug script to test STREME file parsing
"""

import re

def debug_parse_streme_file(streme_file):
    """
    Debug version of STREME parser to see what's happening
    """
    print(f"Parsing file: {streme_file}")
    
    with open(streme_file, 'r') as f:
        content = f.read()
    
    print(f"File content length: {len(content)} characters")
    
    # Look for all MOTIF lines first
    motif_lines = re.findall(r'MOTIF.*', content)
    print(f"Found {len(motif_lines)} MOTIF lines:")
    for i, line in enumerate(motif_lines[:5]):  # Show first 5
        print(f"  {i+1}: {line}")
    
    # Try the new regex pattern
    motif_sections = re.findall(r'MOTIF\s+(\d+)-([ATCGRYSWKMBDHVN]+)\s+(STREME-\d+)(.*?)(?=MOTIF|\Z)', content, re.DOTALL)
    print(f"\nFound {len(motif_sections)} motif sections with new regex")
    
    for i, (motif_number, consensus_seq, streme_id, motif_data) in enumerate(motif_sections[:3]):
        print(f"\nMotif {i+1}:")
        print(f"  Number: {motif_number}")
        print(f"  Consensus: {consensus_seq}")
        print(f"  STREME ID: {streme_id}")
        print(f"  Data length: {len(motif_data)} chars")
        
        # Look for E-value in the data
        evalue_match = re.search(r'E=\s*([0-9.e-]+)', motif_data)
        if evalue_match:
            print(f"  E-value: {evalue_match.group(1)}")
        else:
            print("  E-value: not found")
        
        # Look for nsites
        sites_match = re.search(r'nsites=\s*(\d+)', motif_data)
        if sites_match:
            print(f"  Sites: {sites_match.group(1)}")
        else:
            print("  Sites: not found")

if __name__ == "__main__":
    debug_parse_streme_file("/home/l338m483/scratch/MEME_Test/All_lines/streme_Genes_IM155_DNA_lifted/streme.txt")
