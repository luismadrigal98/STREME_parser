#!/usr/bin/env python3
"""
Validate motif consolidation by analyzing the consolidated_streme_sites.tsv output.
Check if we're inappropriately clustering unrelated sequences together.
"""

import csv
import sys
from collections import defaultdict

def analyze_motif_cluster(motif_id, sites):
    """
    Analyze a single motif cluster to see if sequences are appropriately grouped
    """
    sequences = [site['sequence'] for site in sites]
    streme_patterns = set(site['original_streme_consensus'] for site in sites)
    
    # Calculate sequence diversity metrics
    unique_sequences = set(sequences)
    seq_lengths = [len(seq) for seq in sequences]
    
    # Check base composition variability
    base_compositions = []
    for seq in sequences[:50]:  # Sample first 50 to avoid too much computation
        total_len = len(seq)
        if total_len > 0:
            comp = {
                'A': seq.count('A') / total_len,
                'T': seq.count('T') / total_len,
                'G': seq.count('G') / total_len,
                'C': seq.count('C') / total_len
            }
            base_compositions.append(comp)
    
    # Calculate composition variance (simple measure of diversity)
    if base_compositions:
        avg_A = sum(comp['A'] for comp in base_compositions) / len(base_compositions)
        avg_T = sum(comp['T'] for comp in base_compositions) / len(base_compositions)
        avg_G = sum(comp['G'] for comp in base_compositions) / len(base_compositions)
        avg_C = sum(comp['C'] for comp in base_compositions) / len(base_compositions)
        
        # Calculate variance for each base
        var_A = sum((comp['A'] - avg_A)**2 for comp in base_compositions) / len(base_compositions)
        var_T = sum((comp['T'] - avg_T)**2 for comp in base_compositions) / len(base_compositions)
        var_G = sum((comp['G'] - avg_G)**2 for comp in base_compositions) / len(base_compositions)
        var_C = sum((comp['C'] - avg_C)**2 for comp in base_compositions) / len(base_compositions)
        
        total_variance = var_A + var_T + var_G + var_C
    else:
        total_variance = 0
    
    return {
        'total_sites': len(sites),
        'unique_sequences': len(unique_sequences),
        'sequence_diversity': len(unique_sequences) / len(sites) if sites else 0,
        'streme_patterns_merged': len(streme_patterns),
        'streme_patterns': list(streme_patterns),
        'length_range': (min(seq_lengths), max(seq_lengths)) if seq_lengths else (0, 0),
        'composition_variance': total_variance,
        'sample_sequences': sequences[:5],  # First 5 sequences as examples
        'consensus': sites[0]['motif_consensus'] if sites else '',
        'cluster_size': sites[0]['cluster_size'] if sites else 0
    }

def main():
    if len(sys.argv) != 2:
        print("Usage: python3 validate_consolidation.py consolidated_streme_sites.tsv")
        sys.exit(1)
    
    tsv_file = sys.argv[1]
    
    print("=== MOTIF CONSOLIDATION VALIDATION ===")
    print(f"Analyzing: {tsv_file}")
    print()
    
    # Read the consolidated data
    motif_clusters = defaultdict(list)
    
    try:
        with open(tsv_file, 'r') as f:
            reader = csv.DictReader(f, delimiter='\t')
            for row in reader:
                motif_id = row['consolidated_motif_id']
                motif_clusters[motif_id].append(row)
    except Exception as e:
        print(f"Error reading file: {e}")
        sys.exit(1)
    
    print(f"Found {len(motif_clusters)} consolidated motifs")
    print()
    
    # Analyze each cluster
    suspicious_clusters = []
    
    for motif_id, sites in motif_clusters.items():
        analysis = analyze_motif_cluster(motif_id, sites)
        
        # Flag potentially problematic clusters
        is_suspicious = False
        issues = []
        
        # High sequence diversity might indicate unrelated sequences
        if analysis['sequence_diversity'] > 0.8 and analysis['total_sites'] > 20:
            is_suspicious = True
            issues.append(f"High sequence diversity ({analysis['sequence_diversity']:.2f})")
        
        # Multiple very different STREME patterns merged
        if analysis['streme_patterns_merged'] > 3:
            is_suspicious = True
            issues.append(f"Many STREME patterns merged ({analysis['streme_patterns_merged']})")
        
        # Very high composition variance
        if analysis['composition_variance'] > 0.3:
            is_suspicious = True
            issues.append(f"High composition variance ({analysis['composition_variance']:.3f})")
        
        # Large length range might indicate different motif types
        length_diff = analysis['length_range'][1] - analysis['length_range'][0]
        if length_diff > 5:
            is_suspicious = True
            issues.append(f"Large length range ({analysis['length_range']})")
        
        if is_suspicious:
            suspicious_clusters.append((motif_id, analysis, issues))
    
    # Report results
    print(f"=== VALIDATION RESULTS ===")
    print(f"Total motifs analyzed: {len(motif_clusters)}")
    print(f"Potentially problematic clusters: {len(suspicious_clusters)}")
    print()
    
    if suspicious_clusters:
        print("=== SUSPICIOUS CLUSTERS (May need review) ===")
        
        # Sort by total sites (most common first)
        suspicious_clusters.sort(key=lambda x: x[1]['total_sites'], reverse=True)
        
        for motif_id, analysis, issues in suspicious_clusters[:10]:  # Show top 10
            print(f"\n🚨 {motif_id} ({analysis['total_sites']} sites)")
            print(f"   True consensus: {analysis['consensus']}")
            print(f"   STREME patterns merged: {analysis['streme_patterns']}")
            print(f"   Issues: {', '.join(issues)}")
            print(f"   Sample sequences:")
            for i, seq in enumerate(analysis['sample_sequences'], 1):
                print(f"     {i}. {seq}")
            
            # Recommendation
            if analysis['composition_variance'] > 0.4:
                print("   💡 Recommendation: Very high variance - likely unrelated sequences")
            elif analysis['streme_patterns_merged'] > 5:
                print("   💡 Recommendation: Too many patterns merged - consider stricter threshold")
            else:
                print("   💡 Recommendation: Review manually - may be legitimate variation")
    
    else:
        print("✅ No obviously suspicious clusters found!")
        print("   The consolidation appears to be working appropriately.")
    
    # Summary statistics
    print(f"\n=== SUMMARY STATISTICS ===")
    
    # Overall cluster size distribution
    cluster_sizes = [len(sites) for sites in motif_clusters.values()]
    single_site_clusters = sum(1 for size in cluster_sizes if size == 1)
    small_clusters = sum(1 for size in cluster_sizes if 2 <= size <= 10)
    medium_clusters = sum(1 for size in cluster_sizes if 11 <= size <= 100)
    large_clusters = sum(1 for size in cluster_sizes if size > 100)
    
    print(f"Single-site motifs: {single_site_clusters}")
    print(f"Small clusters (2-10 sites): {small_clusters}")
    print(f"Medium clusters (11-100 sites): {medium_clusters}")
    print(f"Large clusters (>100 sites): {large_clusters}")
    
    # Merged STREME patterns statistics
    streme_merges = [analysis['streme_patterns_merged'] for _, sites in motif_clusters.items() 
                     for analysis in [analyze_motif_cluster(_, sites)]]
    avg_merges = sum(streme_merges) / len(streme_merges) if streme_merges else 0
    max_merges = max(streme_merges) if streme_merges else 0
    
    print(f"Average STREME patterns per cluster: {avg_merges:.1f}")
    print(f"Maximum STREME patterns merged: {max_merges}")
    
    if len(suspicious_clusters) / len(motif_clusters) < 0.05:
        print(f"\n✅ OVERALL ASSESSMENT: Good consolidation quality")
        print(f"   Less than 5% of clusters flagged as suspicious")
    elif len(suspicious_clusters) / len(motif_clusters) < 0.15:
        print(f"\n⚠️  OVERALL ASSESSMENT: Moderate consolidation quality")
        print(f"   Some clusters may need review")
    else:
        print(f"\n🚨 OVERALL ASSESSMENT: Many suspicious clusters")
        print(f"   Consider adjusting consolidation parameters")

if __name__ == "__main__":
    main()
