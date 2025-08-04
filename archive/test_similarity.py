#!/usr/bin/env python3
"""
Test the improved motif similarity function with your examples
"""

import re

# Enhanced similarity function
def calculate_motif_similarity(motif1, motif2):
    """Calculate similarity between two motifs with IUPAC support"""
    seq1 = motif1.get('consensus', '')
    seq2 = motif2.get('consensus', '')
    
    if not seq1 or not seq2:
        return 0.0
    
    # Handle IUPAC nucleotide codes
    iupac_matches = {
        'R': 'AG', 'Y': 'CT', 'S': 'GC', 'W': 'AT', 'K': 'GT', 'M': 'AC',
        'B': 'CGT', 'D': 'AGT', 'H': 'ACT', 'V': 'ACG', 'N': 'ACGT'
    }
    
    def bases_match(base1, base2):
        """Check if two bases can match considering IUPAC codes"""
        if base1 == base2:
            return True
        bases1 = iupac_matches.get(base1, base1)
        bases2 = iupac_matches.get(base2, base2)
        return bool(set(bases1) & set(bases2))
    
    # Try different alignment strategies
    max_similarity = 0.0
    min_len = min(len(seq1), len(seq2))
    max_len = max(len(seq1), len(seq2))
    
    # Direct alignment
    if min_len > 0:
        matches = sum(1 for i in range(min_len) if bases_match(seq1[i], seq2[i]))
        similarity = matches / max_len
        max_similarity = max(max_similarity, similarity)
    
    # Sliding window alignment
    for offset in range(-3, 4):
        matches = 0
        comparisons = 0
        
        for i in range(min_len):
            pos1 = i
            pos2 = i + offset
            
            if 0 <= pos1 < len(seq1) and 0 <= pos2 < len(seq2):
                if bases_match(seq1[pos1], seq2[pos2]):
                    matches += 1
                comparisons += 1
        
        if comparisons >= min_len * 0.7:
            similarity = matches / comparisons
            length_penalty = min_len / max_len
            adjusted_similarity = similarity * length_penalty
            max_similarity = max(max_similarity, adjusted_similarity)
    
    return max_similarity

# Test with your examples
test_motifs = [
    {'consensus': 'AAAAAAAAAAAAAA', 'id': 'IM155_1', 'line': 'IM155'},
    {'consensus': 'AAAAAAAAAAAAAABA', 'id': 'IM444_1', 'line': 'IM444'},
    {'consensus': 'AATTACAACCACCCCY', 'id': 'IM155_2', 'line': 'IM155'},
    {'consensus': 'AATTACAACCACCCC', 'id': 'IM444_3', 'line': 'IM444'},
]

print("Testing motif similarity function:")
print("="*60)

for i in range(len(test_motifs)):
    for j in range(i+1, len(test_motifs)):
        motif1 = test_motifs[i]
        motif2 = test_motifs[j]
        similarity = calculate_motif_similarity(motif1, motif2)
        
        print(f"{motif1['id']}: {motif1['consensus']}")
        print(f"{motif2['id']}: {motif2['consensus']}")
        print(f"Similarity: {similarity:.3f}")
        print(f"Would cluster at 0.6 threshold: {'YES' if similarity >= 0.6 else 'NO'}")
        print(f"Would cluster at 0.8 threshold: {'YES' if similarity >= 0.8 else 'NO'}")
        print("-" * 40)

print("\nExpected results:")
print("- IM155_1 vs IM444_1 (poly-A): HIGH similarity")
print("- IM155_2 vs IM444_3 (AATTACA...): VERY HIGH similarity")
