"""
Author: Mohammad Abbas Naqvi & Avneesh Saravanapavan
Date: 6th March 2025

Worked together on the entire file
"""
import pandas as pd
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


def knn_calc(testing_seq, training_seqs, training_tags, kmer, number_neighbors):
    simu_list = []
    for i in range(len(training_seqs)):
        train_seq = training_seqs[i]
        sim = spectrum(testing_seq, train_seq, kmer)
        simu_list.append((sim, training_tags[i]))

    simu_list.sort(key=lambda x: x[0], reverse=True)
    tippytop_neigh = simu_list[:number_neighbors]

    vote_count = {}
    vote_sum = {}

    for sim, tag in tippytop_neigh:
        if tag not in vote_count:
            vote_count[tag] = 0
            vote_sum[tag] = 0
        vote_count[tag] += 1
        vote_sum[tag] += sim

    max_votes = max(vote_count.values())
    potentials = [label for label, count in vote_count.items() if
                  count == max_votes]
    if len(potentials) == 1:
        return potentials[0]
    else:
        return max(potentials, key=lambda score: vote_sum[score])



if __name__ == "__main__":


    train_length = [10, 30, 100]
    kmer_sizes = [2, 4, 6, 8]
    number_neighbors = [1, 3, 5, 7]

    test_data = parser("./my_project2_copy/test.fasta")
    test_seqs = []
    test_labels = []
    for key, item in test_data.items():
        if "/class=" in key:
            label = key.split("/class=")[1].split()[0].lower()
        else:
            label = "na"
        test_seqs.append(item)
        test_labels.append(label)

    results = []

    for each in train_length:
        print("\n")
        print("Training on size " + str(each) + ":")
        train_exons = parser("./my_project2_copy/train-exons" + str(each) + ".fasta")
        train_introns = parser("./my_project2_copy/train-introns" + str(each) + ".fasta")


        train_seqs = []
        train_tags = []
        for header, seq in train_exons.items():
            train_seqs.append(seq)
            train_tags.append("exon")
        for header, seq in train_introns.items():
            train_seqs.append(seq)
            train_tags.append("intron")

        for kmer in kmer_sizes:
            for num_neighbor in number_neighbors:
                correct = 0
                total = len(test_seqs)
                for i, test_seq in enumerate(test_seqs):
                    pred = knn_calc(test_seq, train_seqs, train_tags, kmer,
                               num_neighbor)
                    if pred == test_labels[i]:
                        correct += 1
                accuracy = correct / total if total > 0 else 0
                results.append((each, kmer, num_neighbor, accuracy))
                print("Neighbors: " + str(num_neighbor) + ", KMER_Count: " + str(kmer) + ", Accuracy: " + str(round(accuracy, 4)))
    df = pd.DataFrame(results, columns=["TrainSize", "KMER", "Neighbors", "Accuracy"])
    df["Accuracy"] = df["Accuracy"].round(4)
    print("\nSummary of Results:")
    print(df)

