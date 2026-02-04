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

def structure_kmer_alignment_text(aln, sts, line_chars=60, out_name='alignment.md', save=True):
        
    rng, aln_len, qry_cov, idnt, aln_sc, aln_sc_norm = sts

    # Formatting options for pandoc md to pdf

    header = '\n'.join(['---',
                        'geometry: margin=1.5cm',
                        'papersize: letter',
                        'sansfont:',
                        'fontsize: 12pt',
                        'urlcolor: blue',
                        'toc:',
                        'toc-depth: 4',
                        '---'])
    
    # Format alignment ranges
    
    rng_table = ['## Alignment ranges',
                 '',
                 '| **query_start** | **query_end** | **subject_start** | **subject_end** |',
                 '| :---: | :---: | :---: | :---: |']
    
    rng_table += [f'| {qs} | {qe} | {ss} | {se} |' for qs,qe,ss,se in rng]
    
    rng_table = '\n'.join(rng_table)
    
    # Format stats
    
    stats_table = ['## Stats',
                   '',
                   '| | |',
                   '| :--- | :---: |']
    
    stats_table += [f'| *alignment length* | {aln_len} |',
                    f'| *query coverage* | {qry_cov:.3f}% |',
                    f'| *identity* | {idnt:.3f}% |',
                    f'| *matches* | {(aln[1,] == "|").sum()} |',
                    f'| *mismatches* | {(aln[1,] == " ").sum()} |',
                    f'| *gaps* | {(aln[1,] == "-").sum()} |',
                    f'| *alignment score* | {aln_sc} |',
                    f'| *normalized alignment score* | {aln_sc_norm:.3f} |']
    
    stats_table = '\n'.join(stats_table)
    
    # Format alignment
    
    pos_max_str_len = len(str(aln.shape[1]))
    
    aln_txt = ['## Alignment',
               '',
               f'query 1-{aln.shape[1] + 1}bp',
               '',
               '```latex']
    
    for i in range(0, aln.shape[1], line_chars):
        
        start_pos = str(i + 1)
        start_pos = ' ' * (pos_max_str_len - len(start_pos)) + start_pos + ' '
        
        end_pos = ' ' + str(min(aln.shape[1], i + line_chars + 1))
        
        block = [start_pos + ''.join(aln[0, i : i + line_chars]) + end_pos,
                 ' ' * (pos_max_str_len + 1) + ''.join(aln[1, i : i + line_chars]),
                 ' ' * (pos_max_str_len + 1) + ''.join(aln[2, i : i + line_chars]),
                 '']
        
        aln_txt.extend(block)
    
    aln_txt.append('```')
    
    aln_txt = '\n'.join(aln_txt)
    
    # Finalize outut text, then print
    
    out_txt = '\n\n'.join([header, rng_table, stats_table, aln_txt])
    
    if save:
        
        with open(out_name, 'w') as out_file:
            
            out_file.write(out_txt)
    
    return out_txt

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

### ------------------MAIN------------------ ###

import numpy as np

from sys import argv

### Import data

query_path = argv[argv.index('--query') + 1]
query_seq = load_fasta(query_path)
query_seq = {name.split('|')[-1] : seq for name,seq in query_seq.items()}

subject_path = argv[argv.index('--subject') + 1]
subject_seq = load_fasta(subject_path)
subject_seq = {name.split('|')[-1] : seq for name,seq in subject_seq.items()}

### Compare query and subject sequences

# Capitalize sequences

query, subject = list(query_seq.values())[0].upper(), list(subject_seq.values())[0].upper()

# Forward alignment

forward_graph, forward_stats = kmer_alignment(query, subject, kmer_size=31, min_anchor_kmers=50, matches_reward=1, mismatches_penalty=-1, gap_penalty=0)

# Reverse alignment

rev_subject = reverse_complementary(subject, 'DNA')

reverse_graph, reverse_stats = kmer_alignment(query, rev_subject, kmer_size=31, min_anchor_kmers=50, matches_reward=1, mismatches_penalty=-1, gap_penalty=0)

if forward_stats[-1] >= reverse_stats[-1] and forward_stats[-1] > 0:
    
    forward_graph_updated, forward_stats_updated = fill_gaps(query, subject, forward_graph, forward_stats)
    
    _ = structure_kmer_alignment_text(forward_graph_updated, forward_stats_updated, line_chars=80, out_name='alignment.md', save=True)

elif reverse_stats[-1] >= forward_stats[-1] and reverse_stats[-1] > 0:
    
    reverse_graph_updated, reverse_stats_updated = fill_gaps(query, rev_subject, reverse_graph, reverse_stats)
    
    _ = structure_kmer_alignment_text(reverse_graph_updated, reverse_stats_updated, line_chars=80, out_name='alignment.md', save=True)

else:

    with open('alignment.md', 'w') as out_file:
            
        out_file.write('')
