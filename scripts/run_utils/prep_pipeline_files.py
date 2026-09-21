#!/usr/bin/env python3

"""
This script preps the necessary files for running the plasmid assembly pipeline
on the UofU CHPC servers

Data is usually provided by the HSC sequencing core in a zip file, containing
1 directory per barcode and a fasta file for reference sequences (optional)
In the fasta file, each sequence name ends with the barcode

After unzipping the file provided by HSC cores, run this script from the root
directory of the pipeline pointing to the unzipped directory with --plasmid_dir
"""

### ---------------------------------------- ###

class split_fasta:
    
    def __init__(self, fasta_path, output_dir='split_fasta'):
        
        self.fasta_path = fasta_path
        
        self.output_dir = output_dir
        
        # Load multifasta
        
        sequences = self.load_fasta(fasta_path)
        
        # Store sequences sizes
        
        self.seq_sizes = {seq_name : len(seq) for seq_name,seq in sequences.items()}
        
        # Create output_directory if missing
        
        makedirs(output_dir, exist_ok=True)
        
        # Save split sequences
        
        self.split_fasta_paths = {}
        
        for seq_name, seq in sequences.items():
            
            if not len(seq):
                
                continue
            
            fasta_out_txt = self.structure_fasta(seq, seq_name, 80)
            
            split_fasta_out_path = f'{output_dir}/{seq_name}.fasta'
            
            self.split_fasta_paths[seq_name] = split_fasta_out_path

            with open(split_fasta_out_path, 'w') as fasta_out_file:
                
                fasta_out_file.write(fasta_out_txt)
    
    ### ------------------------------------ ###
    
    @staticmethod
    def load_fasta(path):
        
        fasta = {}
        
        for chrom in open(path).read().split('>'):
            
            if not len(chrom):
                
                continue
            
            chrom = chrom.split('\n')
            
            chrom_name = chrom[0].split(' ')[0].replace('|', '_')
            
            chrom_seq = ''.join(chrom[1:])
            
            fasta[chrom_name] = chrom_seq
        
        return fasta

    ### ---------------------------------------- ###

    @staticmethod
    def structure_fasta(seq, seq_name='Seq', line_chars=80):
            
        fasta = [f'>{seq_name}']
            
        for i in range(0, len(seq), line_chars):
                
            fasta.append(seq[i : i + line_chars])
            
        fasta = '\n'.join(fasta)
            
        return fasta

### ------------------MAIN------------------ ###

import gzip

from os import chdir, getcwd, listdir, makedirs
from os.path import isdir, exists
from sys import argv

### Output dirs

WORK_DIR = getcwd()

PIPELINE_DIR = WORK_DIR

SPLIT_FASTA_DIR = f'{WORK_DIR}/split_fasta'

SAMPLE_MANIFESTS_DIR = f'{WORK_DIR}/sample_manifests'

REPORTS_DIR = f'{WORK_DIR}/reports'

BASE_SLURM_SCRIPT = """#!/bin/bash
#SBATCH --nodes=1
#SBATCH --account=ucgd-rw
#SBATCH --partition=ucgd-rw
#SBATCH --qos ucgd-prod-rw
#SBATCH -o assembly_[PLASMID_ID]-out-%j
#SBATCH -e assembly_[PLASMID_ID]-err-%j

# Loading modules
module load singularity/4.1.1 openjdk/23.0.1 miniconda3/25.9.1

# Main paths
sample_manifest_path=[MANIFEST_PATH]
main_output_dir=[MAIN_OUTPUT_DIR]

# Running Nextflow
nextflow_executable=/uufs/chpc.utah.edu/common/HIPAA/u1084359/work_dir/nextflow/nextflow

${nextflow_executable} run main.nf \
-profile singularity \
--sample_manifest_path ${sample_manifest_path} \
--main_output_dir $(basename ${main_output_dir})

cp assembly_[PLASMID_ID]-{out,err}-${SLURM_JOB_ID} ${main_output_dir}/

# Generate report
report_dir="report_$(basename ${main_output_dir})"

/bin/bash scripts/run_utils/prep_report.sh ${main_output_dir} ${report_dir}

zip -r ${report_dir}.zip ${report_dir}/*

rm -r ${report_dir}
"""

### Parse plasmid directory

plasmid_dir_path = argv[argv.index('--plasmid_dir') + 1]

plasmid_name = plasmid_dir_path.split('/')[-1]

if not exists(plasmid_dir_path):
    
    print('ERROR: plasmid_dir path invalid.')
    
else:
    
    # Anatomize plasmid_dir path
    
    if WORK_DIR not in plasmid_dir_path:
        
        chdir(plasmid_dir_path)
        
        plasmid_dir_path = getcwd()
        
        chdir(WORK_DIR)
    
    # Find barcode subdirectories
    
    barcode_subdirs = [f for f in listdir(plasmid_dir_path) if f.startswith('barcode') and isdir(f'{plasmid_dir_path}/{f}')]
    
    # Remove barcode_subdirs with too few reads

    reads_files_suffix = '.fastq.gz'

    min_reads = 10
    
    barcodes_to_remove = []
    
    for b in barcode_subdirs:
        
        b_reads = 0
        
        try:
        
            fq_files = [f'{plasmid_dir_path}/{b}/{file}' for file in listdir(f'{plasmid_dir_path}/{b}') if file.endswith(reads_files_suffix)]
            
            for fq in fq_files:
                
                with gzip.open(fq, 'rt') as fq_open:
    
                    b_reads += (sum(1 for line in fq_open) // 4)
            
            if b_reads < min_reads:
                
                barcodes_to_remove.append(b)
    
        except:
            
            # Better to process a bad barcode than to exclude a good one
            continue
    
    barcode_subdirs = [b for b in barcode_subdirs if b not in barcodes_to_remove]
    
    # Find fasta file
    
    fasta_path = [f'{plasmid_dir_path}/{f}' for f in listdir(plasmid_dir_path) if f.endswith('fasta')]
    
    if len(fasta_path):
        
        fasta_path = fasta_path[0]
        
        fasta_info = split_fasta(fasta_path, f'{SPLIT_FASTA_DIR}/{plasmid_name}')
        
        split_fasta_paths = fasta_info.split_fasta_paths

        seq_sizes = fasta_info.seq_sizes
    
    else:
        
        split_fasta_paths = {}
        
        seq_sizes = {}
    
    # Add missing info
    
    for b in barcode_subdirs:
            
        b_fasta_name = [k for k in seq_sizes.keys() if f'_{b}_' in k]
        
        if len(b_fasta_name):
            
            b_fasta_name = b_fasta_name[0]
            
            b_fasta_size = seq_sizes[b_fasta_name]
                
            if b_fasta_size > 0:
                
                # Fasta exists and sequence is not empty
                
                continue
            
            else:
                
                # Sequence is empty, checking for expected size in fasta header
                
                try:
                    
                    seq_sizes[b_fasta_name] = int(b_fasta_name.split('_')[-1])
                
                except:
                    
                    seq_sizes[b_fasta_name] = 10000
            
        else:
                
            split_fasta_paths[b] = 'mock.fasta'
            
            seq_sizes[b] = 10000
    
    # Generate sample manifest
    
    sample_manifest = ['sample_id\tsample_path\tplasmid_size\treference_fasta']
    
    for b in barcode_subdirs:
            
        b_fasta_name = [k for k in seq_sizes.keys() if f'_{b}_' in k]
            
        if len(b_fasta_name):
            
            b_fasta_name = b_fasta_name[0]

            b_fasta_path = split_fasta_paths[b_fasta_name] if b_fasta_name in split_fasta_paths.keys() else ''
            
            b_fasta_size = seq_sizes[b_fasta_name] if b_fasta_name in seq_sizes.keys() else ''
                
            new_entry = '\t'.join([b, f'{plasmid_dir_path}/{b}', str(b_fasta_size), b_fasta_path])
            
        else:
            
            new_entry = '\t'.join([b, f'{plasmid_dir_path}/{b}', '', 'mock.fasta'])
            
        sample_manifest.append(new_entry)
    
    sample_manifest = '\n'.join(sample_manifest)
    
    makedirs(SAMPLE_MANIFESTS_DIR, exist_ok=True)
    
    manifest_out_path = f'{SAMPLE_MANIFESTS_DIR}/{plasmid_name}_manifest.tsv'
    
    with open(manifest_out_path, 'w') as manifest_out:
        
        manifest_out.write(sample_manifest)
    
    # Generate slurm script
    
    slurm_script = BASE_SLURM_SCRIPT
    
    slurm_script = slurm_script.replace('[PLASMID_ID]', plasmid_name)
    
    slurm_script = slurm_script.replace('[MANIFEST_PATH]', manifest_out_path)
    
    slurm_script = slurm_script.replace('[MAIN_OUTPUT_DIR]', f'{WORK_DIR}/assembly_{plasmid_name}')
    
    with open(f'{PIPELINE_DIR}/run_{plasmid_name}_assembly.slurm', 'w') as slurm_script_out:
        
        slurm_script_out.write(slurm_script)
