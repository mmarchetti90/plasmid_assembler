#!/usr/bin/env python3

### ---------------------------------------- ###

def load_fasta(path):
    
    fasta = {}
    
    for chrom in open(path).read().split('>'):
        
        if not len(chrom):
            
            continue
        
        chrom = chrom.split('\n')
        
        chrom_name = chrom[0].split(' ')[0]
        
        chrom_seq = ''.join(chrom[1:])
        
        fasta[chrom_name] = chrom_seq
    
    return fasta

### ---------------------------------------- ###

def structure_fasta(seq, seq_name='Seq', line_chars=80):
        
    fasta = [f'>{seq_name}']
        
    for i in range(0, len(seq), line_chars):
            
        fasta.append(seq[i : i + line_chars])
        
    fasta = '\n'.join(fasta)
        
    return fasta

### ---------------------------------------- ###

def kmer_alignment(query, subject, kmer_size=31, min_anchor_kmers=50, matches_reward=1, mismatches_penalty=-1, gap_penalty=0):
    
    # All uppercase
    
    query, subject = query.upper(), subject.upper()
    
    # Find shared kmers as anchors
    # N.B. using ".index()" when computing match_offsets to make algorithm greedy
    
    query_kmers = [query[i : i + kmer_size] for i in range(0, len(query) - kmer_size + 1, 1)]
    
    subject_kmers = [subject[i : i + kmer_size] for i in range(0, len(subject) - kmer_size + 1, 1)]
    
    match_offsets = np.array([(n, subject_kmers.index(qk), subject_kmers.index(qk) - n)
                              for n,qk in enumerate(query_kmers)
                              if qk in subject_kmers])
    
    # Align
    
    if not len(match_offsets):
        
        return np.array([[], [], []]), ([], 0, 0, 0, 0, 0)
    
    else:
        
        # Align using anchors
        
        alignment = np.repeat('-', 3 * len(query)).reshape((3, len(query)))
        
        alignment[0,] = list(query)
        
        ranges = []
        
        for o in np.unique(match_offsets[:, 2]):
            
            matches_sub = match_offsets[match_offsets[:, 2] == o,].copy()
            
            matches_sub = matches_sub[np.argsort(matches_sub[:, 0])]
            
            if matches_sub.shape[0] < min_anchor_kmers:
                
                continue
            
            q_start, q_end = matches_sub[0, 0], matches_sub[-1, 0]
            
            s_start, s_end = matches_sub[0, 1], matches_sub[-1, 1]
            
            alignment[2, q_start : q_end + kmer_size] = list(subject[s_start : s_end + kmer_size])
        
            ranges.append((q_start, q_end + kmer_size, s_start, s_end + kmer_size))
        
        if not len(ranges):
            
            return np.array([[], [], []]), ([], 0, 0, 0, 0, 0)
        
        else:
            
            # Fill in alignment graph
        
            alignment[1, alignment[0,] == alignment[2,]] = '|'
            
            alignment[1, (alignment[0,] != alignment[2,]) & (alignment[2,] != '-')] = ' '
            
            # Alignment stats
            
            alignment_score = (matches_reward * np.sum(alignment[1,] == '|') +
                               mismatches_penalty * np.sum(alignment[1,] == ' ') +
                               gap_penalty * np.sum(alignment[1,] == '-'))
            
            query_coverage = 100 * np.sum(alignment[1,] != '-') / alignment.shape[1]
            
            alignment_length = np.where(alignment[1,] != '-')[0]
            alignment_length = (alignment_length[-1] - alignment_length[0] + 1) if len(alignment_length) else 0
            
            identity = (100 * np.sum(alignment[1,] == '|') / alignment_length) if alignment_length else 0
            
            norm_alignment_score = alignment_score * (alignment_length / alignment.shape[1])
            
            return alignment, (ranges, alignment_length, query_coverage, identity, alignment_score, norm_alignment_score)

### ---------------------------------------- ###

def reverse_complementary(seq, seq_type='DNA'):
    
    if seq_type == 'DNA':
        
        rev_comp = seq[::-1].lower().replace('a', 'T').replace('t', 'A').replace('g', 'C').replace('c', 'G')
    
    elif seq_type == 'RNA':
        
        rev_comp = seq[::-1].lower().replace('a', 'U').replace('u', 'A').replace('g', 'C').replace('c', 'G')
    
    else:
        
        rev_comp = ''
        
    return rev_comp

### ---------------------------------------- ###

def fill_gaps(qry, sbj, graph, stats, kmer_size=5, min_anchor_kmers=5, matches_reward=1, mismatches_penalty=-1, gap_penalty=0):
    
    rng, aln_len, qry_cov, idnt, aln_sc, aln_sc_norm = stats
    
    # Sort ranges
    
    rng.sort(key=lambda r: r[0])
    
    # Fill gaps between anchors
    
    for r1, r2 in zip(rng, rng[1:]):
        
        if r1[1] < r2[0]:
            
            qry_sub = qry[r1[1] : r2[0]]
            
            sbj_sub = sbj[r1[3] : r2[2]]
    
            sub_graph, _ = kmer_alignment(qry_sub, sbj_sub, kmer_size, min_anchor_kmers)
            
            if sub_graph.shape[1] == r2[0] - r1[1]:
            
                graph[:, r1[1] : r2[0]] = sub_graph
    
    # Alignment stats
    
    alignment_score = (matches_reward * np.sum(graph[1,] == '|') +
                       mismatches_penalty * np.sum(graph[1,] == ' ') +
                       gap_penalty * np.sum(graph[1,] == '-'))
    
    query_coverage = 100 * np.sum(graph[1,] != '-') / graph.shape[1]
    
    alignment_length = np.where(graph[1,] != '-')[0]
    alignment_length = alignment_length[-1] - alignment_length[0] + 1
    
    identity = 100 * np.sum(graph[1,] == '|') / alignment_length
    
    norm_alignment_score = alignment_score * (alignment_length / graph.shape[1])
    
    return graph, (rng, alignment_length, query_coverage, identity, alignment_score, norm_alignment_score)

### ---------------------------------------- ###

def kneedle(vector, sort_vector=True):
    
    """
    Kneedle to find threshold cutoff.
    """
    
    if sort_vector:
        
        vector = np.sort(vector)[::-1]
    
    # Find gradient and intercept
    x0, x1 = 0, len(vector)
    y0, y1 = max(vector), min(vector)
    gradient = (y1 - y0) / (x1 - x0)
    intercept = y0
    
    # Compute difference vector
    difference_vector = [(gradient * x + intercept) - y for x,y in enumerate(vector)]
    
    # Find max of difference_vector and define cutoff
    cutoff_index = difference_vector.index(max(difference_vector))
    cutoff_value = vector[cutoff_index]
    
    return cutoff_index, cutoff_value

### ------------------MAIN------------------ ###

import numpy as np
import pandas as pd

from os import listdir
from sys import argv

### Import fasta contigs

sample_id = argv[argv.index('--sample_id') + 1]

fasta_dir = argv[argv.index('--fasta_dir') + 1]

contigs = [contig.upper() for file in listdir(fasta_dir) if file.endswith('fasta') for contig in load_fasta(f'{fasta_dir}/{file}').values()]

### Compute pairwise alignments

pairwise_scores = {'a' : [],
                   'b' : [],
                   'score' : []}

for query_idx in range(len(contigs) - 1):
    
    query = contigs[query_idx]
    
    for subject_idx in range(query_idx + 1, len(contigs)):
        
        subject = contigs[subject_idx]

        # Forward alignment

        _, stats = kmer_alignment(query, subject, kmer_size=31, min_anchor_kmers=50, matches_reward=1, mismatches_penalty=-1, gap_penalty=0)
        
        forward_score = stats[-1]
        
        # Reverse alignment
        
        rev_subject = reverse_complementary(subject, 'DNA')
        
        _, stats = kmer_alignment(query, rev_subject, kmer_size=31, min_anchor_kmers=50, matches_reward=1, mismatches_penalty=-1, gap_penalty=0)
        
        reverse_score = stats[-1]
        
        # Store scores
        
        pairwise_scores['a'].append(query_idx)
        pairwise_scores['b'].append(subject_idx)
        pairwise_scores['score'].append(max(forward_score, reverse_score))

pairwise_scores = pd.DataFrame(pairwise_scores)

pairwise_scores = pairwise_scores.sort_values(by='score', ascending=False)

pairwise_scores_backup = pairwise_scores.copy()

### Threshold (positive) pairwise alignment scores

pairwise_scores = pairwise_scores.loc[pairwise_scores['score'] > 0,]

_, scores_cutoff = kneedle(pairwise_scores['score'].values, False)

pairwise_scores = pairwise_scores.loc[pairwise_scores['score'] >= scores_cutoff,]

### Cluster based on similarity

# Init clustering with starting pair

assigned_elements, clusters = [], []

while pairwise_scores.shape[0]:
    
    # Create new cluster

    a, b = pairwise_scores.iloc[0, :2].astype(int)
    
    assigned_elements.extend([a, b])
    
    # Update cluster
    
    cluster = [a, b]
    
    update_toggle = True
    
    while update_toggle:
        
        scores_sub = pairwise_scores.loc[((pairwise_scores['a'].isin(cluster)) |
                                          (pairwise_scores['b'].isin(cluster))) &
                                         ((~ pairwise_scores['a'].isin(assigned_elements)) |
                                          (~ pairwise_scores['b'].isin(assigned_elements))),]
        
        if scores_sub.shape[0] > 0:
            
            new_cluster_elements = [element for element in np.unique(scores_sub[['a', 'b']].values.ravel()) if element not in cluster]
            
            assigned_elements.extend(new_cluster_elements)
            
            cluster.extend(new_cluster_elements)
            
            update_toggle = True
        
        else:
            
            update_toggle = False

    clusters.append(cluster)

    # Remove cluster elements from pairwise_scores
    
    pairwise_scores = pairwise_scores.loc[(~ pairwise_scores['a'].isin(assigned_elements)) |
                                          (~ pairwise_scores['b'].isin(assigned_elements)),]

### Discard clusters with less than N elements

N = 3

clusters = [cl for cl in clusters if len(cl) >= N]

### Consensus assembly

assemblies, multiple_sequence_alignements = [], []

for cl_num,cl in enumerate(clusters):
    
    query = contigs[cl[0]]
    
    alignments = [np.array(list(query))]
    
    for subject_idx in cl[1:]:
        
        subject = contigs[subject_idx]
        
        # Forward alignment

        forward_alignment, forward_stats = kmer_alignment(query, subject, kmer_size=31, min_anchor_kmers=50, matches_reward=1, mismatches_penalty=-1, gap_penalty=0)
        
        if forward_alignment.shape[1] > 0:
            
            forward_alignment_updated, forward_stats_updated = fill_gaps(query, subject, forward_alignment, forward_stats)
        
        else:
            
            forward_alignment_updated, forward_stats_updated = forward_alignment, forward_stats
        
        # Reverse alignment
        
        rev_subject = reverse_complementary(subject, 'DNA')
        
        reverse_alignment, reverse_stats = kmer_alignment(query, rev_subject, kmer_size=31, min_anchor_kmers=50, matches_reward=1, mismatches_penalty=-1, gap_penalty=0)
        
        if reverse_alignment.shape[1] > 0:
        
            reverse_alignment_updated, reverse_stats_updated = fill_gaps(query, rev_subject, reverse_alignment, reverse_stats)
        
        else:
            
            reverse_alignment_updated, reverse_stats_updated = reverse_alignment, reverse_stats
        
        # Store data
        
        if forward_stats_updated[-1] > reverse_stats_updated[-1]:
            
            alignments.append(forward_alignment_updated[-1])
        
        else:
            
            alignments.append(reverse_alignment_updated[-1])
        
    alignments = np.stack(alignments)
    
    # Create consensus

    consensus, coverage = [], []

    for base in alignments.T:
        
        base_scores = [(nucleotide, (base == nucleotide).sum()) for nucleotide in 'ATGC']
        
        base_scores.sort(key=lambda bs: bs[1], reverse=True)
        
        nucleotide, score = base_scores[0]
        
        consensus.append(nucleotide)
        
        coverage.append(score)

    consensus = ''.join(consensus)
    
    consensus = structure_fasta(consensus, f'>{sample_id}_CL={cl_num}_N={len(cl)}', 80)

    assemblies.append(consensus)
    
    coverage_flags = ''.join([' ' if cov == len(cl) else '*' for cov in coverage])
    
    # Generate multiple sequence alignment
    
    msa = f'>{sample_id}_CL={cl_num}_N={len(cl)}' + '\n' + '\n'.join([''.join(aln) for aln in alignments]) + '\n' + coverage_flags
    
    multiple_sequence_alignements.append(msa)

### Save to file

with open(f'{sample_id}_custom.consensus.contigs.fasta', 'w') as output:
    
    output.write('\n\n'.join(assemblies))

with open(f'{sample_id}_custom.subassemblies.alignments.txt', 'w') as output:
    
    output.write('\n\n'.join(multiple_sequence_alignements))
