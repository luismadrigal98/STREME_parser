#!/usr/bin/env python3
"""
Comprehensive script to compile all motifs from multiple STREME runs,
identify common vs unique motifs across lines, and create regulatory maps.

This script processes independent STREME results from each line to:
1. Extract all motifs from all lines
2. Compare motifs to identify common vs unique ones
3. Create positional maps for each line
4. Generate regulatory strings for deep learning
"""

import sys
import os
import re
from pathlib import Path
from collections import defaultdict
import glob

def parse_streme_file(streme_file):
    """
    Parse a single STREME output file to extract motifs.
    """
    motifs = []
    
    try:
        with open(streme_file, 'r') as f:
            content = f.read()
        
        # Look for motif sections - STREME format: "MOTIF 1-AAAAAAAAAAAAAA STREME-1"
        motif_sections = re.findall(r'MOTIF\s+(\d+)-([ATCGRYSWKMBDHVN]+)\s+(STREME-\d+)(.*?)(?=MOTIF|\Z)', content, re.DOTALL)
        
        for motif_number, consensus_seq, streme_id, motif_data in motif_sections:
            motif_info = {
                'id': f"{motif_number}-{consensus_seq}",
                'consensus': consensus_seq,  # Extract consensus from motif ID
                'streme_id': streme_id,
                'number': int(motif_number),
                'evalue': '',
                'sites': '',
                'width': len(consensus_seq),
                'pwm': [],
                'raw_sequence': consensus_seq
            }
            
            # Extract E-value (handle scientific notation properly)
            evalue_match = re.search(r'E=\s*([0-9.e+-]+)', motif_data)
            if evalue_match:
                try:
                    motif_info['evalue'] = float(evalue_match.group(1))
                except ValueError:
                    # Try alternate E-value format
                    evalue_match2 = re.search(r'E-value\s*=\s*([0-9.e+-]+)', motif_data)
                    if evalue_match2:
                        motif_info['evalue'] = float(evalue_match2.group(1))
                    else:
                        motif_info['evalue'] = 1.0  # Default value
            
            # Extract sites count from nsites
            sites_match = re.search(r'nsites=\s*(\d+)', motif_data)
            if sites_match:
                motif_info['sites'] = int(sites_match.group(1))
            
            # Extract width from w= parameter
            width_match = re.search(r'w=\s*(\d+)', motif_data)
            if width_match:
                motif_info['width'] = int(width_match.group(1))
            
            # Try to extract PWM data if available
            pwm_section = re.search(r'letter-probability matrix:(.*?)(?=MOTIF|\Z)', motif_data, re.DOTALL)
            if pwm_section:
                pwm_lines = pwm_section.group(1).strip().split('\n')
                for line in pwm_lines:
                    line = line.strip()
                    if line and not line.startswith('letter-probability') and not line.startswith('alength'):
                        try:
                            values = [float(x) for x in line.split()]
                            if len(values) == 4:  # A, C, G, T probabilities
                                motif_info['pwm'].append(values)
                        except ValueError:
                            continue
            
            motifs.append(motif_info)
    
    except Exception as e:
        print(f"Error parsing {streme_file}: {e}")
        return []
    
    return motifs

def calculate_motif_similarity(motif1, motif2):
    """
    Calculate similarity between two motifs based on consensus sequence.
    Returns similarity score (0-1).
    Enhanced for biological motif comparison.
    """
    seq1 = motif1.get('consensus', '')
    seq2 = motif2.get('consensus', '')
    
    if not seq1 or not seq2:
        return 0.0
    
    # Handle IUPAC nucleotide codes for better biological accuracy
    iupac_matches = {
        'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC',
        'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'
    }
    
    def bases_match(base1, base2):
        """Check if two bases can match considering IUPAC codes"""
        if base1 == base2:
            return True
        
        # Get possible bases for each position
        bases1 = iupac_matches.get(base1, base1)
        bases2 = iupac_matches.get(base2, base2)
        
        # Check if there's any overlap
        return bool(set(bases1) & set(bases2))
    
    # Try different alignment strategies
    max_similarity = 0.0
    
    # Strategy 1: Direct alignment
    min_len = min(len(seq1), len(seq2))
    max_len = max(len(seq1), len(seq2))
    
    if min_len > 0:
        matches = sum(1 for i in range(min_len) if bases_match(seq1[i], seq2[i]))
        similarity = matches / max_len  # Penalize length differences
        max_similarity = max(max_similarity, similarity)
    
    # Strategy 2: Sliding window alignment (for offset motifs)
    for offset in range(-3, 4):  # Try small offsets
        matches = 0
        comparisons = 0
        
        for i in range(min_len):
            pos1 = i
            pos2 = i + offset
            
            if 0 <= pos1 < len(seq1) and 0 <= pos2 < len(seq2):
                if bases_match(seq1[pos1], seq2[pos2]):
                    matches += 1
                comparisons += 1
        
        if comparisons >= min_len * 0.7:  # At least 70% overlap
            similarity = matches / comparisons
            # Apply length penalty
            length_penalty = min_len / max_len
            adjusted_similarity = similarity * length_penalty
            max_similarity = max(max_similarity, adjusted_similarity)
    
    # Strategy 3: Check reverse complement
    rev_comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
                'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M', 'M': 'K',
                'B': 'V', 'D': 'H', 'H': 'D', 'V': 'B'}
    
    seq2_rc = ''.join(rev_comp.get(base, base) for base in reversed(seq2))
    
    # Direct RC comparison
    if min_len > 0:
        matches = sum(1 for i in range(min_len) if bases_match(seq1[i], seq2_rc[i]))
        similarity = matches / max_len
        max_similarity = max(max_similarity, similarity)
    
    # Sliding RC comparison
    for offset in range(-3, 4):
        matches = 0
        comparisons = 0
        
        for i in range(min_len):
            pos1 = i
            pos2 = i + offset
            
            if 0 <= pos1 < len(seq1) and 0 <= pos2 < len(seq2_rc):
                if bases_match(seq1[pos1], seq2_rc[pos2]):
                    matches += 1
                comparisons += 1
        
        if comparisons >= min_len * 0.7:
            similarity = matches / comparisons
            length_penalty = min_len / max_len
            adjusted_similarity = similarity * length_penalty
            max_similarity = max(max_similarity, adjusted_similarity)
    
    return max_similarity

def cluster_similar_motifs(all_motifs_by_line, similarity_threshold=0.8):
    """
    Cluster similar motifs across lines and identify common vs unique motifs.
    """
    # Flatten all motifs with line information
    all_motifs = []
    for line_name, motifs in all_motifs_by_line.items():
        for motif in motifs:
            motif['source_line'] = line_name
            all_motifs.append(motif)
    
    # Create clusters of similar motifs
    clusters = []
    assigned = set()
    
    for i, motif1 in enumerate(all_motifs):
        if i in assigned:
            continue
        
        cluster = [i]
        assigned.add(i)
        
        for j, motif2 in enumerate(all_motifs[i+1:], i+1):
            if j in assigned:
                continue
            
            similarity = calculate_motif_similarity(motif1, motif2)
            if similarity >= similarity_threshold:
                cluster.append(j)
                assigned.add(j)
        
        clusters.append(cluster)
    
    # Analyze clusters
    motif_catalog = {
        'common_motifs': [],
        'line_specific_motifs': [],
        'clusters': []
    }
    
    for cluster_idx, cluster in enumerate(clusters):
        cluster_motifs = [all_motifs[i] for i in cluster]
        lines_represented = set(m['source_line'] for m in cluster_motifs)
        
        cluster_info = {
            'cluster_id': f'cluster_{cluster_idx:03d}',
            'motifs': cluster_motifs,
            'lines': list(lines_represented),
            'is_common': len(lines_represented) > 1,
            'representative_motif': cluster_motifs[0],  # Use first as representative
            'consensus_sequences': list(set(m.get('consensus', '') for m in cluster_motifs if m.get('consensus')))
        }
        
        motif_catalog['clusters'].append(cluster_info)
        
        if cluster_info['is_common']:
            motif_catalog['common_motifs'].append(cluster_info)
        else:
            motif_catalog['line_specific_motifs'].append(cluster_info)
    
    return motif_catalog

def create_regulatory_maps(motif_catalog, all_motifs_by_line):
    """
    Create regulatory maps showing motif positions for each line.
    """
    regulatory_maps = {}
    
    for line_name, motifs in all_motifs_by_line.items():
        # Sort motifs by some positional criteria (using E-value as proxy for importance)
        sorted_motifs = sorted(motifs, key=lambda x: x.get('evalue', 1.0))
        
        line_map = {
            'line': line_name,
            'total_motifs': len(motifs),
            'motif_positions': [],
            'regulatory_string': '',
            'motif_types': defaultdict(int)
        }
        
        # Assign each motif to its cluster
        for pos, motif in enumerate(sorted_motifs):
            # Find which cluster this motif belongs to
            cluster_assignment = None
            for cluster in motif_catalog['clusters']:
                for cluster_motif in cluster['motifs']:
                    if (motif.get('id') == cluster_motif.get('id') and 
                        motif.get('consensus') == cluster_motif.get('consensus')):
                        cluster_assignment = cluster['cluster_id']
                        break
                if cluster_assignment:
                    break
            
            motif_entry = {
                'position': pos,
                'motif_id': motif.get('id', f'motif_{pos}'),
                'cluster_id': cluster_assignment or 'unclustered',
                'consensus': motif.get('consensus', ''),
                'evalue': motif.get('evalue', ''),
                'sites': motif.get('sites', ''),
                'is_common': cluster_assignment in [c['cluster_id'] for c in motif_catalog['common_motifs']]
            }
            
            line_map['motif_positions'].append(motif_entry)
            
            # Count motif types
            motif_type = 'common' if motif_entry['is_common'] else 'unique'
            line_map['motif_types'][motif_type] += 1
        
        # Create regulatory string (simplified)
        regulatory_elements = []
        for motif_entry in line_map['motif_positions']:
            # Use cluster ID as regulatory element
            element = motif_entry['cluster_id'][:8]  # Truncate for readability
            regulatory_elements.append(element)
        
        line_map['regulatory_string'] = '-'.join(regulatory_elements)
        regulatory_maps[line_name] = line_map
    
    return regulatory_maps

def write_motif_catalog(motif_catalog, output_file):
    """
    Write comprehensive motif catalog to file.
    """
    with open(output_file, 'w') as f:
        f.write("# Comprehensive Motif Catalog\n")
        f.write("# Generated from all STREME line analyses\n")
        f.write("# " + "="*80 + "\n\n")
        
        # Summary statistics
        f.write("## SUMMARY STATISTICS\n")
        f.write(f"Total motif clusters: {len(motif_catalog['clusters'])}\n")
        f.write(f"Common motifs (present in multiple lines): {len(motif_catalog['common_motifs'])}\n")
        f.write(f"Line-specific motifs: {len(motif_catalog['line_specific_motifs'])}\n\n")
        
        # Common motifs
        f.write("## COMMON MOTIFS (present in multiple lines)\n")
        f.write("ClusterID\tLines\tConsensus\tRepresentativeMotif\n")
        for cluster in motif_catalog['common_motifs']:
            lines_str = ','.join(cluster['lines'])
            consensus_str = ','.join(cluster['consensus_sequences'])
            f.write(f"{cluster['cluster_id']}\t{lines_str}\t{consensus_str}\t{cluster['representative_motif']['id']}\n")
        f.write("\n")
        
        # Line-specific motifs
        f.write("## LINE-SPECIFIC MOTIFS\n")
        f.write("ClusterID\tLine\tConsensus\tMotifID\tE-value\n")
        for cluster in motif_catalog['line_specific_motifs']:
            for motif in cluster['motifs']:
                f.write(f"{cluster['cluster_id']}\t{motif['source_line']}\t{motif.get('consensus', '')}\t{motif['id']}\t{motif.get('evalue', '')}\n")
        f.write("\n")
        
        # Detailed cluster information
        f.write("## DETAILED CLUSTER INFORMATION\n")
        for cluster in motif_catalog['clusters']:
            f.write(f"\n### Cluster: {cluster['cluster_id']}\n")
            f.write(f"Lines represented: {', '.join(cluster['lines'])}\n")
            f.write(f"Number of motifs: {len(cluster['motifs'])}\n")
            f.write(f"Type: {'Common' if cluster['is_common'] else 'Line-specific'}\n")
            f.write("Motifs in cluster:\n")
            for motif in cluster['motifs']:
                f.write(f"  - {motif['id']} ({motif['source_line']}) - {motif.get('consensus', '')} - E-val: {motif.get('evalue', '')}\n")
            f.write("-" * 50 + "\n")

def write_regulatory_maps(regulatory_maps, output_file):
    """
    Write regulatory maps for each line.
    """
    with open(output_file, 'w') as f:
        f.write("# Regulatory Maps by Line\n")
        f.write("# Position and order of motifs in each line (5' to 3')\n")
        f.write("# " + "="*80 + "\n\n")
        
        for line_name, line_map in regulatory_maps.items():
            f.write(f"## LINE: {line_name}\n")
            f.write(f"Total motifs: {line_map['total_motifs']}\n")
            f.write(f"Common motifs: {line_map['motif_types']['common']}\n")
            f.write(f"Unique motifs: {line_map['motif_types']['unique']}\n")
            f.write(f"Regulatory string: {line_map['regulatory_string']}\n\n")
            
            f.write("Motif positions (5' to 3'):\n")
            f.write("Pos\tMotifID\tClusterID\tType\tConsensus\tE-value\n")
            
            for motif in line_map['motif_positions']:
                motif_type = 'Common' if motif['is_common'] else 'Unique'
                f.write(f"{motif['position']}\t{motif['motif_id']}\t{motif['cluster_id']}\t{motif_type}\t{motif['consensus']}\t{motif['evalue']}\n")
            
            f.write("\n" + "="*60 + "\n\n")

def write_regulatory_strings_fasta(regulatory_maps, output_file):
    """
    Write regulatory strings as FASTA for deep learning applications.
    """
    with open(output_file, 'w') as f:
        f.write("# Regulatory Strings for Deep Learning\n# Each sequence represents the motif composition of a line\n")
        
        for line_name, line_map in regulatory_maps.items():
            # Create a more sophisticated regulatory string
            # Using cluster IDs as "amino acids" in regulatory alphabet
            regulatory_sequence = line_map['regulatory_string'].replace('-', '')
            
            f.write(f">{line_name}_regulatory_string\n")
            f.write(f"{regulatory_sequence}\n")

def main():
    if len(sys.argv) < 2:
        print("Usage: python motif_compiler_comprehensive.py <streme_results_directory> [similarity_threshold]")
        print("Example: python motif_compiler_comprehensive.py /path/to/streme_results 0.8")
        sys.exit(1)
    
    results_dir = sys.argv[1]
    similarity_threshold = float(sys.argv[2]) if len(sys.argv) > 2 else 0.8
    
    if not Path(results_dir).exists():
        print(f"Error: Directory {results_dir} does not exist")
        sys.exit(1)
    
    print(f"Processing STREME results from: {results_dir}")
    print(f"Using similarity threshold: {similarity_threshold}")
    
    # Find all STREME result directories
    streme_dirs = glob.glob(os.path.join(results_dir, "*streme*"))
    streme_dirs.extend(glob.glob(os.path.join(results_dir, "*", "streme.txt")))
    
    if not streme_dirs:
        print("No STREME result directories found. Looking for streme.txt files...")
        streme_files = glob.glob(os.path.join(results_dir, "**", "streme.txt"), recursive=True)
        streme_dirs = [os.path.dirname(f) for f in streme_files]
    
    print(f"Found {len(streme_dirs)} STREME result directories")
    
    # Parse all STREME files
    all_motifs_by_line = {}
    
    for streme_dir in streme_dirs:
        line_name = os.path.basename(streme_dir)
        if line_name.startswith('streme_'):
            line_name = line_name[7:]  # Remove 'streme_' prefix
        
        streme_file = os.path.join(streme_dir, 'streme.txt')
        if not os.path.exists(streme_file):
            print(f"Warning: streme.txt not found in {streme_dir}")
            continue
        
        print(f"Processing line: {line_name}")
        motifs = parse_streme_file(streme_file)
        
        if motifs:
            all_motifs_by_line[line_name] = motifs
            print(f"  Found {len(motifs)} motifs")
        else:
            print(f"  No motifs found")
    
    if not all_motifs_by_line:
        print("No motifs found in any files. Exiting.")
        sys.exit(1)
    
    print(f"\nTotal lines processed: {len(all_motifs_by_line)}")
    total_motifs = sum(len(motifs) for motifs in all_motifs_by_line.values())
    print(f"Total motifs found: {total_motifs}")
    
    # Cluster similar motifs
    print("\nClustering similar motifs...")
    motif_catalog = cluster_similar_motifs(all_motifs_by_line, similarity_threshold)
    
    print(f"Created {len(motif_catalog['clusters'])} motif clusters")
    print(f"Common motifs: {len(motif_catalog['common_motifs'])}")
    print(f"Line-specific motifs: {len(motif_catalog['line_specific_motifs'])}")
    
    # Create regulatory maps
    print("\nCreating regulatory maps...")
    regulatory_maps = create_regulatory_maps(motif_catalog, all_motifs_by_line)
    
    # Write outputs
    print("\nWriting output files...")
    
    # Write motif catalog
    write_motif_catalog(motif_catalog, "comprehensive_motif_catalog.txt")
    print("  - comprehensive_motif_catalog.txt")
    
    # Write regulatory maps
    write_regulatory_maps(regulatory_maps, "regulatory_maps_by_line.txt")
    print("  - regulatory_maps_by_line.txt")
    
    # Write regulatory strings for deep learning
    write_regulatory_strings_fasta(regulatory_maps, "regulatory_strings_for_ML.fasta")
    print("  - regulatory_strings_for_ML.fasta")
    
    print("\nAnalysis complete!")
    print("\nFiles generated:")
    print("1. comprehensive_motif_catalog.txt - Complete motif analysis with common vs unique motifs")
    print("2. regulatory_maps_by_line.txt - Positional maps of motifs for each line")
    print("3. regulatory_strings_for_ML.fasta - Regulatory sequences for deep learning")

if __name__ == "__main__":
    main()
