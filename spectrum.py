"""
Author: Mohammad Abbas Naqvi & Avneesh Saravanapavan
Date: 6th March 2025

Worked together on the entire file
"""


def spectrum(seq1, seq2, k):
    kmer_counts1 = {}
    for i in range(len(seq1) - k + 1):
        kmer = seq1[i:i + k]
        kmer_counts1[kmer] = kmer_counts1.get(kmer, 0) + 1
    kmer_counts2 = {}
    for i in range(len(seq2) - k + 1):
        kmer = seq2[i:i + k]
        kmer_counts2[kmer] = kmer_counts2.get(kmer, 0) + 1
    similarity = 0
    for kmer, count in kmer_counts1.items():
        similarity += count * kmer_counts2.get(kmer, 0)

    return similarity