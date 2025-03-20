"""
Author: Mohammad Abbas Naqvi & Avneesh Saravanapavan
Date: 6th March 2025

Worked together on the entire file
"""

Nucleotides = ['A', 'G', 'C', 'T']
Neighbours = {}

def mismatch(sequence1, sequence2, k):
    feature1 = {}
    feature2 = {}

    seq1_len = len(sequence1)
    seq2_len = len(sequence2)

    for i in range(seq1_len - k + 1):
        mykmer = sequence1[i:i+k]
        neighbour = set()
        if (mykmer, 1) in  Neighbours:
            neighbour = Neighbours[(mykmer, 1)]
        else:
            myset = set()
            myset.add(mykmer)
            for pos in range(len(mykmer)):
                for nuc in  Nucleotides:
                    if nuc != mykmer[pos]:
                        temp = mykmer[:pos] + nuc + mykmer[pos+1:]
                        myset.add(temp) # computing new neighbours first sequence
            Neighbours[(mykmer, 1)] = myset
            neighbour = myset
        for dudelivingnexttome in neighbour:
            feature1[dudelivingnexttome] = feature1.get(dudelivingnexttome, 0) + 1
    featurevector1 = feature1


    for i in range(seq2_len - k + 1):
        mykmer2 = sequence2[i:i+k]
        neighbour2 = set()
        if (mykmer2, 1) in  Neighbours:
            neighbour2 = Neighbours[(mykmer2, 1)]
        else:
            myset2 = set()
            myset2.add(mykmer2)
            for pos in range(len(mykmer2)):
                for nuc in  Nucleotides:
                    if nuc != mykmer2[pos]:
                        temp2 = mykmer2[:pos] + nuc + mykmer2[pos+1:]
                        myset2.add(temp2) # computing new neighbours second sequence
            Neighbours[(mykmer2, 1)] = myset2
            neighbour2 = myset2
        for dudelivingnexttome in neighbour2:
            feature2[dudelivingnexttome] = feature2.get(dudelivingnexttome, 0) + 1
    featurevector2 = feature2

    mykernelval = 0
    for key,value in featurevector1.items():
        if key in featurevector2.keys():
            mykernelval += value * featurevector2[key]

    return mykernelval