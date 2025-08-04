#!/usr/bin/env python3
"""
Multi-line motif analysis for Mimulus guttatus lines.
Compares motifs across lines, identifies common/unique motifs,
and creates regulatory strings for deep learning applications.
"""

import sys
import os
import re
from pathlib import Path
from collections import defaultdict
import json

def extract_line_id(directory_name):
    """Extract line ID from directory name like 'streme_Genes_IM1034_DNA_lifted'"""
    match = re.search(r'IM(\d+)', directory_name)
    return match.group(1) if match else directory_name

def parse_streme_output(streme_file):
    """
    Parse STREME output file and extract motif information.
    """
    motifs = []
    
    if not os.path.exists(streme_file):
        print(f"Warning: STREME file not found: {streme_file}")
        return motifs
    
    with open(streme_file, 'r') as f:
        content = f.read()
    
    # Parse motifs using regex patterns
    motif_blocks = re.findall(r'MOTIF\s+(\S+).*?(?=MOTIF|\Z)', content, re.DOTALL)
    
    for i, block in enumerate(motif_blocks):
        lines = block.split('\n')
        motif_info = {
            'id': f'motif_{i+1}',
            'consensus': '',
            'evalue': 'N/A',
            'sites': 'N/A',
            'width': 'N/A',
            'pwm': [],
            'sequence_logos': []
        }
        
        # Extract motif name from first line
        first_line = lines[0].strip()
        parts = first_line.split()
        if len(parts) > 0:
            motif_info['id'] = parts[0]
        
        # Look for E-value, sites, and other info in the block
        for line in lines:
            # E-value
            evalue_match = re.search(r'E-value[:\s=]+([0-9.e-]+)', line, re.IGNORECASE)
            if evalue_match:
                motif_info['evalue'] = evalue_match.group(1)
            
            # Sites
            sites_match = re.search(r'(\d+)\s+sites', line, re.IGNORECASE)
            if sites_match:
                motif_info['sites'] = sites_match.group(1)
            
            # Width
            width_match = re.search(r'width[:\s=]+(\d+)', line, re.IGNORECASE)
            if width_match:
                motif_info['width'] = width_match.group(1)
        
        motifs.append(motif_info)
    
    return motifs

def calculate_motif_similarity(motif1, motif2):
    """
    Simple motif similarity calculation based on available information.
    This is a placeholder - you might want to implement more sophisticated
    similarity measures based on PWMs or consensus sequences.
    """
    # Compare by width first
    try:
        width1 = int(motif1.get('width', 0))
        width2 = int(motif2.get('width', 0))
        
        if abs(width1 - width2) > 3:  # Width difference threshold
            return 0.0
        
        # Compare E-values (lower is better)
        try:
            eval1 = float(motif1.get('evalue', 1.0))
            eval2 = float(motif2.get('evalue', 1.0))
            
            # Simple similarity based on E-value ranges
            if eval1 < 0.001 and eval2 < 0.001:
                return 0.8  # Both highly significant
            elif eval1 < 0.01 and eval2 < 0.01:
                return 0.6  # Both significant
            elif eval1 < 0.05 and eval2 < 0.05:
                return 0.4  # Both moderately significant
            else:
                return 0.2  # Low similarity
        except:
            return 0.3  # Default similarity if E-values not comparable
            
    except:
        return 0.1  # Very low similarity if comparison fails

def find_similar_motifs_across_lines(all_motifs, similarity_threshold=0.5):
    """
    Find motifs that are similar across different lines.
    """
    motif_clusters = []
    used_motifs = set()
    
    lines = list(all_motifs.keys())
    
    for line1 in lines:
        for motif1 in all_motifs[line1]:
            motif1_key = f"{line1}_{motif1['id']}"
            
            if motif1_key in used_motifs:
                continue
                
            # Start a new cluster
            cluster = {
                'cluster_id': f"cluster_{len(motif_clusters) + 1}",
                'motifs': [(line1, motif1)],
                'lines': [line1],
                'consensus_info': {
                    'avg_width': int(motif1.get('width', 0)) if motif1.get('width', '0').isdigit() else 0,
                    'best_evalue': motif1.get('evalue', 'N/A'),
                    'total_sites': int(motif1.get('sites', 0)) if motif1.get('sites', '0').isdigit() else 0
                }
            }
            used_motifs.add(motif1_key)
            
            # Find similar motifs in other lines
            for line2 in lines:
                if line2 == line1:
                    continue
                    
                for motif2 in all_motifs[line2]:
                    motif2_key = f"{line2}_{motif2['id']}"
                    
                    if motif2_key in used_motifs:
                        continue
                    
                    similarity = calculate_motif_similarity(motif1, motif2)
                    
                    if similarity >= similarity_threshold:
                        cluster['motifs'].append((line2, motif2))
                        cluster['lines'].append(line2)
                        used_motifs.add(motif2_key)
                        
                        # Update cluster consensus info
                        if motif2.get('width', '0').isdigit():
                            cluster['consensus_info']['avg_width'] = (
                                cluster['consensus_info']['avg_width'] + int(motif2['width'])
                            ) // 2
                        
                        if motif2.get('sites', '0').isdigit():
                            cluster['consensus_info']['total_sites'] += int(motif2['sites'])
            
            motif_clusters.append(cluster)
    
    return motif_clusters

def create_regulatory_map(all_motifs, motif_clusters):
    """
    Create regulatory maps for each line showing motif presence and order.
    """
    regulatory_maps = {}
    
    for line_id, motifs in all_motifs.items():
        regulatory_maps[line_id] = {
            'line_id': line_id,
            'total_motifs': len(motifs),
            'motifs': [],
            'motif_order': [],
            'common_motifs': [],
            'unique_motifs': []
        }
        
        # Assign each motif to clusters or mark as unique
        for motif in motifs:
            motif_entry = {
                'motif_id': motif['id'],
                'evalue': motif.get('evalue', 'N/A'),
                'sites': motif.get('sites', 'N/A'),
                'width': motif.get('width', 'N/A'),
                'cluster_id': None,
                'is_common': False
            }
            
            # Find which cluster this motif belongs to
            for cluster in motif_clusters:
                for cluster_line, cluster_motif in cluster['motifs']:
                    if cluster_line == line_id and cluster_motif['id'] == motif['id']:
                        motif_entry['cluster_id'] = cluster['cluster_id']
                        motif_entry['is_common'] = len(cluster['lines']) > 1
                        break
                if motif_entry['cluster_id']:
                    break
            
            regulatory_maps[line_id]['motifs'].append(motif_entry)
            
            if motif_entry['is_common']:
                regulatory_maps[line_id]['common_motifs'].append(motif_entry)
            else:
                regulatory_maps[line_id]['unique_motifs'].append(motif_entry)
        
        # Create motif order (5' to 3')
        regulatory_maps[line_id]['motif_order'] = [m['motif_id'] for m in regulatory_maps[line_id]['motifs']]
    
    return regulatory_maps

def create_regulatory_strings(regulatory_maps, motif_clusters):
    """
    Create regulatory strings for deep learning applications.
    """
    regulatory_strings = {}
    
    # Create a mapping from cluster IDs to simplified identifiers
    cluster_map = {}
    for i, cluster in enumerate(motif_clusters):
        cluster_map[cluster['cluster_id']] = f"M{i+1:02d}"
    
    for line_id, reg_map in regulatory_maps.items():
        # Method 1: Cluster-based string
        cluster_string = ""
        for motif in reg_map['motifs']:
            if motif['cluster_id']:
                cluster_string += cluster_map[motif['cluster_id']] + "-"
            else:
                cluster_string += "UNQ-"  # Unique motif
        
        # Method 2: Significance-based string
        significance_string = ""
        for motif in reg_map['motifs']:
            try:
                evalue = float(motif['evalue']) if motif['evalue'] != 'N/A' else 1.0
                if evalue < 0.001:
                    significance_string += "H"  # High significance
                elif evalue < 0.01:
                    significance_string += "M"  # Medium significance
                elif evalue < 0.05:
                    significance_string += "L"  # Low significance
                else:
                    significance_string += "N"  # Not significant
            except:
                significance_string += "N"
        
        # Method 3: Combined regulatory string
        combined_string = ""
        for motif in reg_map['motifs']:
            motif_code = ""
            
            # Add cluster info
            if motif['cluster_id']:
                motif_code += cluster_map[motif['cluster_id']]
            else:
                motif_code += "UNQ"
            
            # Add significance
            try:
                evalue = float(motif['evalue']) if motif['evalue'] != 'N/A' else 1.0
                if evalue < 0.001:
                    motif_code += "H"
                elif evalue < 0.01:
                    motif_code += "M"
                else:
                    motif_code += "L"
            except:
                motif_code += "L"
            
            combined_string += motif_code + "_"
        
        regulatory_strings[line_id] = {
            'cluster_string': cluster_string.rstrip('-'),
            'significance_string': significance_string,
            'combined_string': combined_string.rstrip('_'),
            'motif_count': len(reg_map['motifs']),
            'common_motif_count': len(reg_map['common_motifs']),
            'unique_motif_count': len(reg_map['unique_motifs'])
        }
    
    return regulatory_strings

def write_comprehensive_output(all_motifs, motif_clusters, regulatory_maps, regulatory_strings, output_dir):
    """
    Write comprehensive analysis results to multiple output files.
    """
    os.makedirs(output_dir, exist_ok=True)
    
    # 1. Motif catalog per line
    with open(f"{output_dir}/motifs_per_line.txt", 'w') as f:
        f.write("# Motifs Per Line Analysis\n")
        f.write("# Format: Line_ID\tMotif_ID\tE-value\tSites\tWidth\tCluster_ID\tIs_Common\n")
        f.write("# " + "="*80 + "\n\n")
        
        for line_id, reg_map in regulatory_maps.items():
            f.write(f"## Line: {line_id}\n")
            f.write(f"Total motifs: {reg_map['total_motifs']}\n")
            f.write(f"Common motifs: {len(reg_map['common_motifs'])}\n")
            f.write(f"Unique motifs: {len(reg_map['unique_motifs'])}\n\n")
            
            for motif in reg_map['motifs']:
                f.write(f"{line_id}\t{motif['motif_id']}\t{motif['evalue']}\t")
                f.write(f"{motif['sites']}\t{motif['width']}\t")
                f.write(f"{motif['cluster_id'] or 'UNIQUE'}\t{motif['is_common']}\n")
            f.write("\n")
    
    # 2. Motif positions and relative order
    with open(f"{output_dir}/motif_positions_and_order.txt", 'w') as f:
        f.write("# Motif Positions and Relative Order (5' to 3')\n")
        f.write("# Format: Line_ID\tPosition\tMotif_ID\tCluster_ID\tE-value\n")
        f.write("# " + "="*80 + "\n\n")
        
        for line_id, reg_map in regulatory_maps.items():
            f.write(f"## Line: {line_id}\n")
            f.write("Motif order (5' to 3'): " + " -> ".join(reg_map['motif_order']) + "\n\n")
            
            for pos, motif in enumerate(reg_map['motifs'], 1):
                f.write(f"{line_id}\t{pos}\t{motif['motif_id']}\t")
                f.write(f"{motif['cluster_id'] or 'UNIQUE'}\t{motif['evalue']}\n")
            f.write("\n")
    
    # 3. Regulatory strings for deep learning
    with open(f"{output_dir}/regulatory_strings_for_ML.txt", 'w') as f:
        f.write("# Regulatory Strings for Deep Learning\n")
        f.write("# These strings encode motif presence, order, and significance\n")
        f.write("# " + "="*80 + "\n\n")
        
        f.write("## Format Explanation:\n")
        f.write("# Cluster String: Motif clusters separated by '-'\n")
        f.write("# Significance String: H=High(E<0.001), M=Medium(E<0.01), L=Low(E<0.05), N=Not significant\n")
        f.write("# Combined String: Cluster+Significance separated by '_'\n\n")
        
        for line_id, strings in regulatory_strings.items():
            f.write(f"Line: {line_id}\n")
            f.write(f"  Cluster String:      {strings['cluster_string']}\n")
            f.write(f"  Significance String: {strings['significance_string']}\n")
            f.write(f"  Combined String:     {strings['combined_string']}\n")
            f.write(f"  Motif Stats: Total={strings['motif_count']}, ")
            f.write(f"Common={strings['common_motif_count']}, Unique={strings['unique_motif_count']}\n\n")
    
    # 4. Motif clusters summary
    with open(f"{output_dir}/motif_clusters.txt", 'w') as f:
        f.write("# Motif Clusters Across Lines\n")
        f.write("# Shows which motifs are common vs. unique between lines\n")
        f.write("# " + "="*80 + "\n\n")
        
        for cluster in motif_clusters:
            f.write(f"Cluster: {cluster['cluster_id']}\n")
            f.write(f"  Lines with this motif: {', '.join(cluster['lines'])}\n")
            f.write(f"  Number of lines: {len(cluster['lines'])}\n")
            f.write(f"  Average width: {cluster['consensus_info']['avg_width']}\n")
            f.write(f"  Total sites: {cluster['consensus_info']['total_sites']}\n")
            f.write("  Motifs in cluster:\n")
            
            for line, motif in cluster['motifs']:
                f.write(f"    {line}: {motif['id']} (E={motif['evalue']}, sites={motif['sites']})\n")
            f.write("\n")
    
    # 5. FASTA file with regulatory strings for deep learning
    with open(f"{output_dir}/regulatory_sequences_for_ML.fasta", 'w') as f:
        for line_id, strings in regulatory_strings.items():
            f.write(f">Line_{line_id}_regulatory_string\n")
            f.write(f"{strings['combined_string']}\n")
    
    # 6. JSON summary for programmatic access
    summary_data = {
        'analysis_summary': {
            'total_lines': len(all_motifs),
            'total_clusters': len(motif_clusters),
            'lines_analyzed': list(all_motifs.keys())
        },
        'regulatory_maps': regulatory_maps,
        'motif_clusters': motif_clusters,
        'regulatory_strings': regulatory_strings
    }
    
    with open(f"{output_dir}/analysis_summary.json", 'w') as f:
        json.dump(summary_data, f, indent=2)

def main():
    if len(sys.argv) != 3:
        print("Usage: python analyze_multi_line_motifs.py <streme_directories_path> <output_directory>")
        print("\nExample: python analyze_multi_line_motifs.py /home/l338m483/scratch/MEME_Test/All_lines ./motif_analysis_results")
        sys.exit(1)
    
    streme_dir_path = sys.argv[1]
    output_dir = sys.argv[2]
    
    if not os.path.exists(streme_dir_path):
        print(f"Error: Directory {streme_dir_path} does not exist")
        sys.exit(1)
    
    print(f"Analyzing STREME outputs in: {streme_dir_path}")
    print(f"Output will be saved to: {output_dir}")
    
    # Find all STREME directories
    streme_dirs = [d for d in os.listdir(streme_dir_path) 
                   if d.startswith('streme_') and os.path.isdir(os.path.join(streme_dir_path, d))]
    
    if not streme_dirs:
        print("No STREME directories found!")
        sys.exit(1)
    
    print(f"Found {len(streme_dirs)} STREME directories:")
    for d in streme_dirs:
        print(f"  {d}")
    
    # Parse motifs from all lines
    all_motifs = {}
    
    for streme_dir in streme_dirs:
        line_id = extract_line_id(streme_dir)
        streme_file = os.path.join(streme_dir_path, streme_dir, 'streme.txt')
        
        print(f"\nProcessing line {line_id} from {streme_file}")
        motifs = parse_streme_output(streme_file)
        all_motifs[line_id] = motifs
        print(f"  Found {len(motifs)} motifs")
    
    if not any(all_motifs.values()):
        print("No motifs found in any STREME output!")
        sys.exit(1)
    
    # Find similar motifs across lines
    print(f"\nFinding similar motifs across lines...")
    motif_clusters = find_similar_motifs_across_lines(all_motifs)
    print(f"Found {len(motif_clusters)} motif clusters")
    
    # Create regulatory maps
    print("Creating regulatory maps...")
    regulatory_maps = create_regulatory_map(all_motifs, motif_clusters)
    
    # Create regulatory strings
    print("Creating regulatory strings for deep learning...")
    regulatory_strings = create_regulatory_strings(regulatory_maps, motif_clusters)
    
    # Write comprehensive output
    print(f"Writing results to {output_dir}...")
    write_comprehensive_output(all_motifs, motif_clusters, regulatory_maps, regulatory_strings, output_dir)
    
    print("\nAnalysis complete! Output files generated:")
    print(f"  {output_dir}/motifs_per_line.txt")
    print(f"  {output_dir}/motif_positions_and_order.txt")
    print(f"  {output_dir}/regulatory_strings_for_ML.txt")
    print(f"  {output_dir}/motif_clusters.txt")
    print(f"  {output_dir}/regulatory_sequences_for_ML.fasta")
    print(f"  {output_dir}/analysis_summary.json")

if __name__ == "__main__":
    main()
