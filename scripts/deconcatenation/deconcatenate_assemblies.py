#!/usr/bin/env python3

### ---------------------------------------- ###

def parse_args():
    
    # Load assembly
    
    assembly_path = argv[argv.index('--assembly') + 1]
    
    assembly = load_fasta(assembly_path)
    
    # Expected assembly size
    
    expected_size = int(argv[argv.index('--expected_size') + 1])
    
    return assembly, expected_size

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

def resolve_concatamers(contig_name, contig, expected_size, kmer_size=1000):
    
    print(f'## Deconcatenating {contig_name}')
    
    # All uppercase
    
    contig = contig.upper()
    
    # Skip if contig is short enough already
    
    toggle = True if len(contig) >= (expected_size * 1.1) else False
    
    if not toggle:
        
        print('Skipping: contig is already close to desired size')
    
    # Detect concatamers

    deconcatenation_pass = 0
    
    min_proof = (len(contig) - expected_size - kmer_size) // 4
    
    while toggle:
        
        deconcatenation_pass += 1
        
        print(f'### Pass {deconcatenation_pass}')
        
        # Deconstruct in kmers
        
        kmers = extract_kmers(contig, kmer_size)
        
        # Find possible concatamers positions
        
        possible_concatamers = {} # Dict with concatamer start as index, and list of rightmost end and proof count as values
        
        for kmer_seq in kmers:
            
            matches_pos = [match.start() for match in re.finditer(kmer_seq, contig)]
            
            base_pos = matches_pos[0]
            
            matches_pos = [mp - base_pos for mp in matches_pos if (mp - base_pos) > 0]
        
            for mp in matches_pos:
                
                if mp in possible_concatamers.keys():
                
                    possible_concatamers[mp][0] = max(possible_concatamers[mp][0], mp + base_pos + kmer_size)
                    
                    possible_concatamers[mp][1] += 1
                
                else:
                    
                    possible_concatamers[mp] = [mp + base_pos + kmer_size, 1]
        
        # Sort and select most likely concatamer start
        
        possible_concatamers = [(start, end, proof) for start,(end,proof) in possible_concatamers.items() if proof >= min_proof]

        possible_concatamers.sort(key=lambda c: c[2], reverse=True)
        
        # Clipping
        
        original_size = len(contig)
        
        if len(possible_concatamers):
            
            concatamer_start, concatamer_end, proof = possible_concatamers[0]
            
            # Clip
            # If the concatamer does not extend from concatamer_start to original_size, then
            # whatever is between concatamer_end and original_size is added to the beginning of the
            # plasmid
            
            contig = contig[concatamer_end:] + contig[:concatamer_start]
            
            # Try another round
            
            toggle = True
            
            print(f'Found concatamers. Original contig size was {original_size}bp, now {len(contig)}bp')
        
        else:
            
            contig = contig[:original_size]
            
            # Stop looking for concatamers
            
            toggle = False
            
            print("Didn't find any more concatamers")
    
    print('')

    return contig

### ---------------------------------------- ###

def extract_kmers(seq, k=21):
    
    kmers = list(set([seq[i : i + k] for i in range(0, len(seq) - k, 1)]))
    
    return kmers

### ---------------------------------------- ###

def structure_fasta(seq, seq_name='Seq', line_chars=80):
        
        fasta = [f'>{seq_name}']
        
        for i in range(0, len(seq), line_chars):
            
            fasta.append(seq[i : i + line_chars])
        
        fasta = '\n'.join(fasta)
        
        return fasta

### ------------------MAIN------------------ ###

import re

from sys import argv

### Parse args and load data

assembly, expected_size = parse_args()

### Find concatamers

deconcatenated_assembly = {}

for a_name, a_seq in assembly.items():
    
    deconcatenated_assembly[a_name] = resolve_concatamers(a_name, a_seq, expected_size, 1000)

### Export fasta

structured_fasta = [structure_fasta(da_seq, da_name) for (da_name,da_seq) in deconcatenated_assembly.items()]

structured_fasta = '\n\n'.join(structured_fasta)

with open('deconcatenated.fasta', 'w') as fasta_out:
    
    fasta_out.write(structured_fasta)
