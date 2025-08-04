#!/usr/bin/env python3

import sys
import os

# Add the path to find motif_compiler_comprehensive.py
sys.path.append('/mnt/1692B2EF92B2D28B/Ongoing_projects/MEME_related/')

from motif_compiler_comprehensive import calculate_motif_similarity

def test_problematic_motifs():
    """Test the two motifs that shouldn't be clustered together"""
    
    motif1 = {
        'consensus': 'AAAAAAAATTA',
        'sites': 15,
        'evalue': 1.4e-02
    }
    
    motif2 = {
        'consensus': 'CTYTCTCTCTCTH',
        'sites': 17,
        'evalue': 1.5e-02
    }
    
    print("Testing problematic motif pair:")
    print(f"Motif 1: {motif1['consensus']} (sites: {motif1['sites']}, evalue: {motif1['evalue']})")
    print(f"Motif 2: {motif2['consensus']} (sites: {motif2['sites']}, evalue: {motif2['evalue']})")
    
    similarity = calculate_motif_similarity(motif1, motif2)
    print(f"Calculated similarity: {similarity:.4f}")
    print(f"Would cluster with threshold 0.8? {'YES' if similarity >= 0.8 else 'NO'}")
    print(f"Would cluster with threshold 0.6? {'YES' if similarity >= 0.6 else 'NO'}")
    print(f"Would cluster with threshold 0.4? {'YES' if similarity >= 0.4 else 'NO'}")
    
    # Let's also test some clearly similar motifs
    print("\n" + "="*50)
    print("Testing clearly similar motifs:")
    
    similar_motif1 = {
        'consensus': 'AAAAAAA',
        'sites': 20,
        'evalue': 1e-03
    }
    
    similar_motif2 = {
        'consensus': 'AAAAAAT',
        'sites': 18,
        'evalue': 2e-03
    }
    
    print(f"Motif A: {similar_motif1['consensus']}")
    print(f"Motif B: {similar_motif2['consensus']}")
    
    similarity2 = calculate_motif_similarity(similar_motif1, similar_motif2)
    print(f"Calculated similarity: {similarity2:.4f}")
    print(f"Would cluster with threshold 0.8? {'YES' if similarity2 >= 0.8 else 'NO'}")

if __name__ == "__main__":
    test_problematic_motifs()
