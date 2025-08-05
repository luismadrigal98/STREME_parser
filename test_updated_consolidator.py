#!/usr/bin/env python3
"""
Test the updated consolidator with overlap merging
"""

# Test just the key functions
import sys
sys.path.append('cli_tools')

# Import the functions we need
from streme_sites_consolidator import merge_overlapping_motifs

# Test data with overlapping sites
test_sites = [
    {
        'consolidated_motif_id': 'MOTIF_001',
        'original_motif_id': '1-AAAAAAAAAAAAAAA',
        'motif_consensus': 'AAAAAAAAAAAAAAA',
        'line': 'IM1034',
        'gene_id': 'MiIM7v11000019m.g',
        'start_pos': 786,
        'end_pos': 800,
        'strand': '-',
        'score': 10.06,
        'sequence': 'AAAATAAAAAAAAAG',
        'representative_sequence': 'AAAAAAAAAAAAAAA',
        'cluster_size': 6,
        'length': 15
    },
    {
        'consolidated_motif_id': 'MOTIF_001',
        'original_motif_id': '1-AAAAAAAAAAAAAAA',
        'motif_consensus': 'AAAAAAAAAAAAAAA',
        'line': 'IM1034',
        'gene_id': 'MiIM7v11000019m.g',
        'start_pos': 787,
        'end_pos': 801,
        'strand': '-',
        'score': 13.84,
        'sequence': 'AAAAATAAAAAAAAA',
        'representative_sequence': 'AAAAAAAAAAAAAAA',
        'cluster_size': 6,
        'length': 15
    },
    {
        'consolidated_motif_id': 'MOTIF_002',
        'original_motif_id': '2-CGCGCGCG',
        'motif_consensus': 'CGCGCGCG',
        'line': 'IM1034',
        'gene_id': 'MiIM7v11000019m.g',
        'start_pos': 1000,
        'end_pos': 1007,
        'strand': '+',
        'score': 8.5,
        'sequence': 'CGCGCGCG',
        'representative_sequence': 'CGCGCGCG',
        'cluster_size': 2,
        'length': 8
    }
]

print("=== Testing Overlap Merging ===")
print(f"Original sites: {len(test_sites)}")
for site in test_sites:
    print(f"  {site['consolidated_motif_id']} at {site['start_pos']}-{site['end_pos']}")

merged = merge_overlapping_motifs(test_sites, overlap_threshold=0.5)

print(f"\nMerged sites: {len(merged)}")
for site in merged:
    merged_count = site.get('merged_count', 1)
    print(f"  {site['consolidated_motif_id']} at {site['start_pos']}-{site['end_pos']}, merged: {merged_count}")

print(f"\nReduction: {len(test_sites)} -> {len(merged)} sites")
