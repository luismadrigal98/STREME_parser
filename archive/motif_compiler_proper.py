#!/usr/bin/env python3
"""
Proper motif compiler using established sequence alignment methods.
"""

import os
import re
import sys
from pathlib import Path
from collections import defaultdict

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

def parse_motif_from_line(line):
    """Parse motif information from STREME output line"""
    # Match the motif ID pattern
    motif_match = re.search(r'(\d+-[A-Z]+)', line)
    if not motif_match:
        return None
    
    motif_id = motif_match.group(1)
    consensus = motif_id.split('-')[1]  # Extract consensus from ID
    
    # Parse E-value (handles scientific notation)
    evalue_match = re.search(r'(\d+\.?\d*e?[+-]?\d*)', line)
    evalue = evalue_match.group(1) if evalue_match else 'N/A'
    
    # Parse sites
    sites_match = re.search(r'(\d+)\s+sites', line)
    sites = sites_match.group(1) if sites_match else 'N/A'
    
    return {
        'id': motif_id,
        'consensus': consensus,
        'evalue': evalue,
        'sites': sites,
        'width': str(len(consensus))
    }

def read_streme_motifs(streme_file):
    """Read motifs from STREME output file"""
    motifs = []
    
    try:
        with open(streme_file, 'r') as f:
            content = f.read()
            
        # Find all MOTIF lines in STREME format: "MOTIF 1-CONSENSUS STREME-1"
        motif_lines = re.findall(r'MOTIF\s+(\d+)-([A-Z]+)\s+STREME-\d+', content)
        
        for motif_match in motif_lines:
            motif_num = motif_match[0]
            consensus = motif_match[1]
            motif_id = f"{motif_num}-{consensus}"
            
            # Try to find the corresponding nsites and E-value for this motif
            # Look for the pattern after the motif declaration
            pattern = rf'MOTIF\s+{re.escape(motif_id)}\s+STREME-\d+.*?nsites=\s*(\d+)\s+E=\s*([\d\.e\-\+]+)'
            match = re.search(pattern, content, re.DOTALL)
            
            if match:
                sites = match.group(1)
                evalue = match.group(2)
            else:
                sites = 'N/A'
                evalue = 'N/A'
            
            motifs.append({
                'id': motif_id,
                'consensus': consensus,
                'evalue': evalue,
                'sites': sites,
                'width': str(len(consensus))
            })
                    
    except Exception as e:
        print(f"Warning: Could not read {streme_file}: {e}")
    
    return motifs

def process_all_lines_sequential(base_dir):
    """
    Process lines sequentially using your suggested logic:
    1. Start with line 1 motifs as unique
    2. For each subsequent line, check if motifs already exist
    3. Only add truly new motifs
    """
    # Find all line directories
    line_dirs = []
    for item in os.listdir(base_dir):
        item_path = os.path.join(base_dir, item)
        if os.path.isdir(item_path) and item.startswith('streme_Genes_'):
            line_dirs.append(item)
    
    line_dirs.sort()
    print(f"Found {len(line_dirs)} line directories: {line_dirs}")
    
    # Track unique motifs
    unique_motifs = []  # List of unique motif dictionaries
    motif_line_occurrence = []  # List of sets tracking which lines have each motif
    line_motif_map = defaultdict(list)  # Which motifs each line has
    
    for line_idx, line_dir in enumerate(line_dirs):
        print(f"\nProcessing line {line_idx + 1}/{len(line_dirs)}: {line_dir}")
        
        # Find STREME output file
        streme_file = os.path.join(base_dir, line_dir, 'streme.txt')
        if not os.path.exists(streme_file):
            print(f"  Warning: No streme.txt found in {line_dir}")
            continue
        
        # Read motifs from this line
        line_motifs = read_streme_motifs(streme_file)
        print(f"  Found {len(line_motifs)} motifs from STREME")
        
        # Process each motif from this line
        for motif in line_motifs:
            print(f"    Checking motif: {motif['consensus']}")
            
            # Check if this motif is similar to any existing unique motif
            found_match = False
            
            for i, existing_motif in enumerate(unique_motifs):
                if sequences_are_similar(motif['consensus'], existing_motif['consensus']):
                    # Found a match - add this line to the occurrence set
                    motif_line_occurrence[i].add(line_dir)
                    line_motif_map[line_dir].append(existing_motif['id'])
                    found_match = True
                    print(f"      -> Matches existing motif {existing_motif['consensus']}")
                    break
            
            if not found_match:
                # New unique motif
                unique_motifs.append(motif)
                motif_line_occurrence.append({line_dir})
                line_motif_map[line_dir].append(motif['id'])
                print(f"      -> NEW unique motif added")
    
    return unique_motifs, motif_line_occurrence, line_motif_map, line_dirs

def write_consolidated_output(unique_motifs, motif_line_occurrence, line_motif_map):
    """Write all output files"""
    
    # 1. Consolidated catalog
    with open('consolidated_motif_catalog_proper.txt', 'w') as f:
        f.write("# Consolidated Motif Catalog (Proper Algorithm)\n")
        f.write("# MotifID\tConsensus\tE-value\tSites\tWidth\tOccurrences\tLines\n")
        
        for i, motif in enumerate(unique_motifs):
            lines = sorted(list(motif_line_occurrence[i]))
            f.write(f"{motif['id']}\t{motif['consensus']}\t{motif['evalue']}\t"
                   f"{motif['sites']}\t{motif['width']}\t{len(lines)}\t{','.join(lines)}\n")
    
    # 2. Regulatory maps by line
    motif_lookup = {motif['id']: motif for motif in unique_motifs}
    
    with open('regulatory_maps_proper.txt', 'w') as f:
        f.write("# Regulatory Maps by Line (Proper Algorithm)\n")
        f.write("# Line\tMotifCount\tMotifs\n")
        
        for line, motif_ids in sorted(line_motif_map.items()):
            motif_consensuses = [motif_lookup[mid]['consensus'] for mid in motif_ids if mid in motif_lookup]
            f.write(f"{line}\t{len(motif_consensuses)}\t{','.join(motif_consensuses)}\n")
    
    # 3. ML strings
    with open('regulatory_strings_proper.fasta', 'w') as f:
        for line, motif_ids in sorted(line_motif_map.items()):
            motif_consensuses = [motif_lookup[mid]['consensus'] for mid in motif_ids if mid in motif_lookup]
            regulatory_string = '|'.join(motif_consensuses)
            f.write(f">{line}\n{regulatory_string}\n")

def main():
    if len(sys.argv) != 2:
        print("Usage: python motif_compiler_proper.py <base_directory>")
        sys.exit(1)
    
    base_dir = sys.argv[1]
    
    if not os.path.exists(base_dir):
        print(f"Error: Directory {base_dir} does not exist")
        sys.exit(1)
    
    print(f"Processing motifs from: {base_dir}")
    print("Using proper sequence comparison with IUPAC handling")
    
    # Process all lines sequentially
    unique_motifs, motif_line_occurrence, line_motif_map, line_dirs = process_all_lines_sequential(base_dir)
    
    print(f"\n=== FINAL RESULTS ===")
    print(f"Total unique motifs found: {len(unique_motifs)}")
    print(f"Lines processed: {len(line_dirs)}")
    
    # Categorize motifs
    common_motifs = []
    line_specific_motifs = []
    
    for i, motif in enumerate(unique_motifs):
        occurrence_count = len(motif_line_occurrence[i])
        if occurrence_count > 1:
            common_motifs.append((motif, motif_line_occurrence[i]))
        else:
            line_specific_motifs.append((motif, motif_line_occurrence[i]))
    
    print(f"Common motifs (>1 line): {len(common_motifs)}")
    print(f"Line-specific motifs: {len(line_specific_motifs)}")
    
    # Write outputs
    write_consolidated_output(unique_motifs, motif_line_occurrence, line_motif_map)
    
    print("\nOutput files generated:")
    print("- consolidated_motif_catalog_proper.txt")
    print("- regulatory_maps_proper.txt") 
    print("- regulatory_strings_proper.fasta")

if __name__ == "__main__":
    main()
