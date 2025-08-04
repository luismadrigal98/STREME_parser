#!/usr/bin/env python3
"""
Test the proper similarity function with problematic motifs.
"""

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

def sequences_are_similar(seq1, seq2, threshold=0.8):
    """
    Check if two sequences are similar enough to be the same motif.
    """
    # Quick filters
    if abs(len(seq1) - len(seq2)) > 3:  # Too different in length
        return False
    
    if len(seq1) < 3 or len(seq2) < 3:  # Too short to compare reliably
        return False
    
    # Use IUPAC-aware comparison
    similarity = expand_iupac_for_comparison(seq1, seq2)
    
    return similarity >= threshold

# Test the problematic motifs
print("=== TESTING PROBLEMATIC MOTIFS ===")
motif1 = "AAAAAAAATTA"
motif2 = "CTYTCTCTCTCTH"

print(f"Motif 1: {motif1} (length: {len(motif1)})")
print(f"Motif 2: {motif2} (length: {len(motif2)})")
print(f"Length difference: {abs(len(motif1) - len(motif2))}")

similarity = expand_iupac_for_comparison(motif1, motif2)
print(f"Similarity score: {similarity:.3f}")
print(f"Should cluster (threshold 0.8)? {sequences_are_similar(motif1, motif2)}")
print("Expected: NO - completely different motifs")

print("\n=== TESTING ACTUALLY SIMILAR MOTIFS ===")

# Test similar AT-rich motifs
similar1 = "AAAAAAAATTA"
similar2 = "AAAAAAATTA"  # One less A
similar3 = "AAAAAAAAAA"  # All A's

print(f"\nMotif A: {similar1}")
print(f"Motif B: {similar2}")
similarity_ab = expand_iupac_for_comparison(similar1, similar2)
print(f"Similarity: {similarity_ab:.3f}")
print(f"Should cluster? {sequences_are_similar(similar1, similar2)}")

print(f"\nMotif A: {similar1}")
print(f"Motif C: {similar3}")
similarity_ac = expand_iupac_for_comparison(similar1, similar3)
print(f"Similarity: {similarity_ac:.3f}")
print(f"Should cluster? {sequences_are_similar(similar1, similar3)}")

print("\n=== TESTING IUPAC CODES ===")
iupac1 = "CTYTCTCTCTCTH"  # Y=C/T, H=A/C/T
iupac2 = "CTCTCTCTCTCTA"  # Concrete version

print(f"IUPAC motif: {iupac1}")
print(f"Concrete:    {iupac2}")
similarity_iupac = expand_iupac_for_comparison(iupac1, iupac2)
print(f"Similarity: {similarity_iupac:.3f}")
print(f"Should cluster? {sequences_are_similar(iupac1, iupac2)}")

print("\n=== TESTING VERY DIFFERENT MOTIFS ===")
diff1 = "AAAAAAA"  # All A's
diff2 = "GGGGGGG"  # All G's

print(f"Motif 1: {diff1}")
print(f"Motif 2: {diff2}")
similarity_diff = expand_iupac_for_comparison(diff1, diff2)
print(f"Similarity: {similarity_diff:.3f}")
print(f"Should cluster? {sequences_are_similar(diff1, diff2)}")
print("Expected: NO - no matching bases")
