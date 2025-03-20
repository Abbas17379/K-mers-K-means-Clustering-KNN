"""
Author: Mohammad Abbas Naqvi & Avneesh Saravanapavan
Date: 6th March 2025

Worked together on the entire file
"""

import random
import math
Nucleotides = ['A', 'G', 'C', 'T']
Neighbours = {}
mylarge = 1.7976931348623157e+308
def parser(fasta_file):
    with open(fasta_file, 'r') as fasta_file:
        seq_dict = {}
        sequence = []
        header = None
        for line in fasta_file:
            line = line.strip()
            if line[0] == ">":
                sequence = []
                header = line[1:]
            else:
                sequence.append(line)
                seq_dict[header] = ''.join(sequence)
    return seq_dict

def spectrum(seq1, seq2, k):
    kmer_counts1 = {}
    for i in range(len(seq1) - k + 1):
        kmer = seq1[i:i+k]
        kmer_counts1[kmer] = kmer_counts1.get(kmer, 0) + 1
    kmer_counts2 = {}
    for i in range(len(seq2) - k + 1):
        kmer = seq2[i:i+k]
        kmer_counts2[kmer] = kmer_counts2.get(kmer, 0) + 1
    similarity = 0
    for kmer, count in kmer_counts1.items():
        similarity += count * kmer_counts2.get(kmer, 0)

    return similarity

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

def k_means(sequences, k, specOrmis, cluster_count):
    sequencelistlen = len(sequences)
    kernelMatrix = [[float(0)] * sequencelistlen for y in range(sequencelistlen)]
    for a in range(sequencelistlen):
        for b in range(a, sequencelistlen):
            if specOrmis == "spec":
                temp = spectrum(sequences[a], sequences[b], k)
            else:
                temp = mismatch(sequences[a], sequences[b], k)
            
            kernelMatrix[a][b] = temp
            kernelMatrix[b][a] = temp

    randomAssignment = [random.randint(0, cluster_count - 1) for x in range(sequencelistlen)]
    for x in range(100):
        my_cluster_dict ={}
        for i in range(cluster_count):
            my_cluster_dict[i] = []
        for i in range(len(randomAssignment)):
            mycluster = randomAssignment[i]
            my_cluster_dict[mycluster].append(i)
        
        recomputed_assignment = []
        for y in range(sequencelistlen):
            idealCluster = None
            idealDist = mylarge
            for z in range(cluster_count):
                #Base case
                if (len(my_cluster_dict[z]) < 1):
                    tempdist = mylarge
                else:
                    kernelsum = float(0)
                    for j in my_cluster_dict[z]:
                        kernelsum += kernelMatrix[y][j]
                    
                    clustersum = float(0)
                    for j in my_cluster_dict[z]:
                        for l in my_cluster_dict[z]:
                            clustersum += kernelMatrix[j][l]

                    tempdist = math.sqrt(kernelMatrix[y][y] - (2.0/len(my_cluster_dict[z]))*kernelsum + (1.0/(len(my_cluster_dict[z])**2))*clustersum)
                    if tempdist < idealDist:
                        idealDist = tempdist
                        idealCluster = z
            recomputed_assignment.append(idealCluster)
        if recomputed_assignment == randomAssignment:
            break
        else:
            randomAssignment = recomputed_assignment
    return recomputed_assignment


if __name__ == "__main__":
    fasta_file = "my_project2_copy/kmeans.fasta"
    seq_dict = parser(fasta_file)
    sequences = []
    classes = []
    for header, sequence in seq_dict.items():
        sequences.append(sequence)
        if "/class=" in header:
            splitHeader = header.split("/class=")  
            afterClass = splitHeader[1] 
            firstWord = afterClass.split()[0]  
            cls = firstWord.strip() 
        else:
            cls = "na"
        classes.append(cls)
    possibleClasses = ["intergenic", "intron", "exon"]
    for kerntype in ["spec", "mismatch"]:
        if kerntype == "spec":
            name = "SPECTRUM"
        else:
            name = "MISMATCH"
        for kmercount in [2, 6]:
            for clustercount in [2, 3, 5]:
                print("\n" + name + " KERNEL (KMER=" + str(kmercount) + ", CLUSTERS=" + str(clustercount) + "):")
                assignments = k_means(sequences, kmercount, kerntype, clustercount)

                clusterTotals = {}  
                for i in range(clustercount):  
                    clusterTotals[i] = 0  

                clusterStats = {}  

                for i in range(clustercount):  
                    clusterStats[i] = {}  

                    for Class in possibleClasses: 
                        clusterStats[i][Class] = 0  


                for i, clusterType in enumerate(assignments):
                    clusterTotals[clusterType] += 1
                    classLabel = classes[i]
                    if classLabel in possibleClasses:
                        clusterStats[clusterType][classLabel] += 1

                for c in range(0, clustercount):
                    print("CLUSTER " + str(c + 1) + ":")
                    total = clusterTotals[c]
                    for cls in possibleClasses:
                        count = clusterStats[c - 1][cls]
                        prop = count / total if total > 0 else 0
                        print(cls + " = " + str(round(prop, 2)) + " (" + str(count) + ")")



            
