#!/usr/bin/env python3
"""
Advanced consolidation script for pair-wise STREME results analysis.
Processes individual line results and identifies similar motifs between pairs.
"""

import os
import sys
import re
from pathlib import Path
from collections import defaultdict
import subprocess

def parse_streme_file(streme_file):
    """
    Parse a STREME output file and extract motif information.
    """
    motifs = []
    
    if not os.path.exists(streme_file):
        return motifs
    
    with open(streme_file, 'r') as f:
        content = f.read()
    
    # Find all motifs using regex
    motif_pattern = r'MOTIF\s+(\S+).*?(?=MOTIF|\Z)'
    motif_matches = re.findall(motif_pattern, content, re.DOTALL)
    
    for i, motif_id in enumerate(motif_matches, 1):
        # Extract more detailed information for each motif
        motif_section = re.search(rf'MOTIF\s+{re.escape(motif_id)}.*?(?=MOTIF|\Z)', content, re.DOTALL)
        
        if motif_section:
            section_text = motif_section.group(0)
            
            # Extract E-value
            evalue_match = re.search(r'E-value[=:\s]+([0-9.e-]+)', section_text)
            evalue = evalue_match.group(1) if evalue_match else 'N/A'
            
            # Extract sites count
            sites_match = re.search(r'(\d+)\s+sites', section_text)
            sites = sites_match.group(1) if sites_match else 'N/A'
            
            # Extract width
            width_match = re.search(r'w=\s*(\d+)', section_text)
            width = width_match.group(1) if width_match else 'N/A'
            
            motifs.append({
                'id': motif_id,
                'index': i,
                'evalue': evalue,
                'sites': sites,
                'width': width,
                'full_section': section_text[:200] + '...' if len(section_text) > 200 else section_text
            })
    
    return motifs

def analyze_line_results(work_dir):
    """
    Analyze all individual line STREME results.
    """
    streme_dirs = sorted([d for d in os.listdir(work_dir) if d.startswith('streme_')])
    
    line_results = {}
    
    for i, streme_dir in enumerate(streme_dirs, 1):
        streme_file = os.path.join(work_dir, streme_dir, 'streme.txt')
        
        print(f"Analyzing Line {i}: {streme_dir}")
        
        motifs = parse_streme_file(streme_file)
        
        line_results[i] = {
            'directory': streme_dir,
            'streme_file': streme_file,
            'motifs': motifs,
            'motif_count': len(motifs)
        }
        
        print(f"  Found {len(motifs)} motifs")
    
    return line_results

def create_pair_analysis(line_results, output_dir):
    """
    Create pair-wise analysis of lines.
    """
    pairs_dir = os.path.join(output_dir, 'detailed_pair_analysis')
    os.makedirs(pairs_dir, exist_ok=True)
    
    line_numbers = sorted(line_results.keys())
    pair_analyses = []
    
    # Process lines in pairs
    for i in range(0, len(line_numbers), 2):
        line1 = line_numbers[i]
        line2 = line_numbers[i + 1] if i + 1 < len(line_numbers) else None
        
        pair_num = (i // 2) + 1
        
        print(f"Creating Pair {pair_num} analysis: Line {line1}" + (f" and Line {line2}" if line2 else " (single)"))
        
        pair_file = os.path.join(pairs_dir, f'pair_{pair_num}_detailed_analysis.txt')
        
        with open(pair_file, 'w') as f:
            f.write(f"# Detailed Pair Analysis: Pair {pair_num}\n")
            f.write(f"# Lines: {line1}" + (f" and {line2}" if line2 else " (single)") + "\n")
            f.write(f"# Generated: {os.path.basename(__file__)}\n")
            f.write("# " + "="*60 + "\n\n")
            
            # Analyze Line 1
            f.write(f"## LINE {line1} DETAILED ANALYSIS\n")
            f.write(f"Source directory: {line_results[line1]['directory']}\n")
            f.write(f"Total motifs: {line_results[line1]['motif_count']}\n\n")
            
            if line_results[line1]['motifs']:
                f.write("Motifs found:\n")
                for motif in line_results[line1]['motifs']:
                    f.write(f"  {motif['index']}. {motif['id']}\n")
                    f.write(f"     E-value: {motif['evalue']}\n")
                    f.write(f"     Sites: {motif['sites']}\n")
                    f.write(f"     Width: {motif['width']}\n\n")
            else:
                f.write("No motifs found.\n\n")
            
            # Analyze Line 2 if exists
            if line2:
                f.write(f"## LINE {line2} DETAILED ANALYSIS\n")
                f.write(f"Source directory: {line_results[line2]['directory']}\n")
                f.write(f"Total motifs: {line_results[line2]['motif_count']}\n\n")
                
                if line_results[line2]['motifs']:
                    f.write("Motifs found:\n")
                    for motif in line_results[line2]['motifs']:
                        f.write(f"  {motif['index']}. {motif['id']}\n")
                        f.write(f"     E-value: {motif['evalue']}\n")
                        f.write(f"     Sites: {motif['sites']}\n")
                        f.write(f"     Width: {motif['width']}\n\n")
                else:
                    f.write("No motifs found.\n\n")
            
            # Comparison section
            f.write("## PAIR COMPARISON\n")
            if line2:
                motifs1 = line_results[line1]['motifs']
                motifs2 = line_results[line2]['motifs']
                
                f.write(f"Line {line1}: {len(motifs1)} motifs\n")
                f.write(f"Line {line2}: {len(motifs2)} motifs\n\n")
                
                # Simple comparison based on motif characteristics
                if motifs1 and motifs2:
                    f.write("Potential similarities (based on width and E-value):\n")
                    
                    similarities = []
                    for m1 in motifs1:
                        for m2 in motifs2:
                            try:
                                # Compare E-values and widths
                                eval1 = float(m1['evalue']) if m1['evalue'] != 'N/A' else 1.0
                                eval2 = float(m2['evalue']) if m2['evalue'] != 'N/A' else 1.0
                                width1 = int(m1['width']) if m1['width'] != 'N/A' else 0
                                width2 = int(m2['width']) if m2['width'] != 'N/A' else 0
                                
                                # Simple similarity criteria
                                if (abs(width1 - width2) <= 2 and 
                                    abs(eval1 - eval2) / max(eval1, eval2, 1e-10) < 0.5):
                                    similarities.append((m1, m2))
                            except:
                                continue
                    
                    if similarities:
                        for m1, m2 in similarities:
                            f.write(f"  {m1['id']} (Line {line1}) ~ {m2['id']} (Line {line2})\n")
                            f.write(f"    Widths: {m1['width']} vs {m2['width']}\n")
                            f.write(f"    E-values: {m1['evalue']} vs {m2['evalue']}\n\n")
                    else:
                        f.write("  No obvious similarities detected.\n\n")
                
                f.write("Recommendations:\n")
                f.write("1. Use TOMTOM to compare motifs between these lines\n")
                f.write("2. Consider merging sequences if motifs are highly similar\n")
                f.write("3. Keep separate if motifs are distinct\n\n")
            else:
                f.write("Single line in this pair - no comparison possible.\n\n")
        
        pair_analyses.append({
            'pair_num': pair_num,
            'line1': line1,
            'line2': line2,
            'file': pair_file
        })
    
    return pair_analyses

def create_master_report(line_results, pair_analyses, output_dir):
    """
    Create a master consolidation report.
    """
    master_file = os.path.join(output_dir, 'master_detailed_report.txt')
    
    with open(master_file, 'w') as f:
        f.write("# Master STREME Consolidation Report\n")
        f.write("# Detailed Analysis of Individual Line Results\n")
        f.write("# " + "="*60 + "\n\n")
        
        # Overall statistics
        total_lines = len(line_results)
        total_motifs = sum(result['motif_count'] for result in line_results.values())
        total_pairs = len(pair_analyses)
        
        f.write("## OVERALL STATISTICS\n")
        f.write(f"Total lines analyzed: {total_lines}\n")
        f.write(f"Total motifs found: {total_motifs}\n")
        f.write(f"Total pairs created: {total_pairs}\n")
        f.write(f"Average motifs per line: {total_motifs/total_lines:.1f}\n\n")
        
        # Per-line summary
        f.write("## PER-LINE SUMMARY\n")
        for line_num in sorted(line_results.keys()):
            result = line_results[line_num]
            f.write(f"Line {line_num}: {result['motif_count']} motifs ({result['directory']})\n")
        f.write("\n")
        
        # Pair summary
        f.write("## PAIR ANALYSIS SUMMARY\n")
        for pair in pair_analyses:
            f.write(f"Pair {pair['pair_num']}: Line {pair['line1']}")
            if pair['line2']:
                f.write(f" & Line {pair['line2']}")
                motifs1 = line_results[pair['line1']]['motif_count']
                motifs2 = line_results[pair['line2']]['motif_count']
                f.write(f" ({motifs1} + {motifs2} motifs)")
            else:
                f.write(" (single)")
                motifs1 = line_results[pair['line1']]['motif_count']
                f.write(f" ({motifs1} motifs)")
            f.write("\n")
        
        f.write("\n## RECOMMENDATIONS\n")
        f.write("1. Review individual pair analyses for detailed comparisons\n")
        f.write("2. Use TOMTOM for precise motif similarity assessment\n")
        f.write("3. Consider consolidating lines with highly similar motifs\n")
        f.write("4. Maintain separate analyses for distinct motif patterns\n\n")
        
        f.write("## FILES GENERATED\n")
        f.write(f"Detailed pair analyses: {len(pair_analyses)} files in detailed_pair_analysis/\n")
        f.write("Master report: master_detailed_report.txt\n")
    
    return master_file

def main():
    if len(sys.argv) != 3:
        print("Usage: python advanced_streme_consolidation.py <work_directory> <output_directory>")
        print("Example: python advanced_streme_consolidation.py /path/to/All_lines ./consolidated_results")
        sys.exit(1)
    
    work_dir = sys.argv[1]
    output_dir = sys.argv[2]
    
    if not os.path.exists(work_dir):
        print(f"Error: Work directory {work_dir} does not exist")
        sys.exit(1)
    
    os.makedirs(output_dir, exist_ok=True)
    
    print("=== Advanced STREME Results Consolidation ===")
    print(f"Work directory: {work_dir}")
    print(f"Output directory: {output_dir}")
    print()
    
    # Analyze individual line results
    print("Step 1: Analyzing individual line results...")
    line_results = analyze_line_results(work_dir)
    print(f"Analyzed {len(line_results)} lines\n")
    
    # Create pair analyses
    print("Step 2: Creating pair-wise analyses...")
    pair_analyses = create_pair_analysis(line_results, output_dir)
    print(f"Created {len(pair_analyses)} pair analyses\n")
    
    # Create master report
    print("Step 3: Creating master report...")
    master_file = create_master_report(line_results, pair_analyses, output_dir)
    print(f"Master report: {master_file}\n")
    
    print("=== Consolidation Complete ===")
    print(f"Results saved to: {output_dir}")
    print("\nNext steps:")
    print("1. Review the master report")
    print("2. Examine individual pair analyses")
    print("3. Run TOMTOM for detailed motif comparisons")
    print("4. Decide which lines to consolidate based on motif similarity")

if __name__ == "__main__":
    main()
