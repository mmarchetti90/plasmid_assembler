#!/bin/bash

### CLI args

input_dir=${1}

output_dir=${2:-assembly_results}

mkdir -p ${output_dir}

### Update permission to avoid issues with singularity

chmod -R 775 ${input_dir}

### Functions

run_pandoc () {

	# N.B. The entrypoint is lost when converting to sif

	in_file=${1}
	out_file=$(basename "${in_file}" | sed "s/.md/.html/g")

	singularity exec \
	-B "$(pwd):/data" \
	/uufs/chpc.utah.edu/common/HIPAA/u1084359/work_dir/singularity_images/pandoc/pandoc.sif \
	pandoc \
	${in_file} -o ${out_file} \
	--syntax-highlighting tango

}

summarize_run () {

	results_dir=${1}

	mkdir run_summaries_tmp

	singularity exec \
	-B "$(pwd):/home" \
	-B ${results_dir} \
	/uufs/chpc.utah.edu/common/HIPAA/u1084359/work_dir/singularity_images/python/python.sif \
	python scripts/run_utils/summarize_individual_barcodes_run.py --results_dir ${results_dir}

	mv barcode*_summary.tsv run_summaries_tmp/

}

### Summarize run

#mkdir run_summaries_tmp

#python scripts/run_utils/summarize_individual_barcodes_run.py --results_dir ${input_dir}

#mv barcode*_summary.tsv run_summaries_tmp/

summarize_run ${input_dir}

### Structure data

for barcode_dir in ${input_dir}/barcode* 
do

	barcode_id=$(basename ${barcode_dir})

	# Reference-based assembly

	mkdir -p ${output_dir}/${barcode_id}/alignment_based_assembly

	cp ${barcode_dir}/alignment_based_assembly/* ${output_dir}/${barcode_id}/alignment_based_assembly/
	mkdir -p ${output_dir}/${barcode_id}/alignment_based_assembly/minimap2_alignment
	cp ${barcode_dir}/minimap2_alignment/*.alignment.log ${output_dir}/${barcode_id}/alignment_based_assembly/minimap2_alignment/
	mkdir -p ${output_dir}/${barcode_id}/alignment_based_assembly/ab1
	cp ${barcode_dir}/ab1/*reference_based_assembly*.ab1 ${output_dir}/${barcode_id}/alignment_based_assembly/ab1/

	# De-novo assembly using all reads

	mkdir -p ${output_dir}/${barcode_id}/de_novo_assembly_all_reads/assembly_annotations
	mkdir -p ${output_dir}/${barcode_id}/de_novo_assembly_all_reads/reference_alignment

	cp ${barcode_dir}/consensus_assembly/*_allreads.polished.contigs.fasta ${output_dir}/${barcode_id}/de_novo_assembly_all_reads/
	cp ${barcode_dir}/annotated_assembly/*_allreads_pLann.{csv,html,gbk} ${output_dir}/${barcode_id}/de_novo_assembly_all_reads/assembly_annotations/
	cp ${barcode_dir}/blast_alignment/*_allreads_blastn.txt ${output_dir}/${barcode_id}/de_novo_assembly_all_reads/reference_alignment/
	run_pandoc ${barcode_dir}/kmer_alignment/${barcode_id}_allreads_alignment.md
	mv ${barcode_id}_allreads_alignment.html ${output_dir}/${barcode_id}/de_novo_assembly_all_reads/reference_alignment/
	mkdir -p ${output_dir}/${barcode_id}/de_novo_assembly_all_reads/ab1
	cp ${barcode_dir}/ab1/*all_reads_assembly*.ab1 ${output_dir}/${barcode_id}/de_novo_assembly_all_reads/ab1/

	# De-novo assembly using reads subsets

	mkdir -p ${output_dir}/${barcode_id}/de_novo_assembly_reads_subsets/intermediate_assemblies
	mkdir -p ${output_dir}/${barcode_id}/de_novo_assembly_reads_subsets/final_assemblies
	mkdir -p ${output_dir}/${barcode_id}/de_novo_assembly_reads_subsets/assembly_annotations
	mkdir -p ${output_dir}/${barcode_id}/de_novo_assembly_reads_subsets/reference_alignment
	mkdir -p ${output_dir}/${barcode_id}/de_novo_assembly_reads_subsets/ab1

	cp ${barcode_dir}/deconcatenated_assembly/*_sample_*_flye.contigs.deconcatenated.fasta ${output_dir}/${barcode_id}/de_novo_assembly_reads_subsets/intermediate_assemblies/
	cp ${barcode_dir}/consensus_assembly/*_cluster_*_trycycler.consensus.polished.contigs.fasta ${output_dir}/${barcode_id}/de_novo_assembly_reads_subsets/final_assemblies/
	cp ${barcode_dir}/consensus_assembly/*_custom.consensus.polished.contigs.fasta ${output_dir}/${barcode_id}/de_novo_assembly_reads_subsets/final_assemblies/
	cp ${barcode_dir}/annotated_assembly/*_{trycycler,custom}.consensus_pLann.{csv,html,gbk} ${output_dir}/${barcode_id}/de_novo_assembly_reads_subsets/assembly_annotations/
	cp ${barcode_dir}/blast_alignment/*_{trycycler,custom}.consensus_blastn.txt ${output_dir}/${barcode_id}/de_novo_assembly_reads_subsets/reference_alignment/
	for kmer_alignment in ${barcode_dir}/kmer_alignment/*_{trycycler,custom}.consensus_alignment.md
	do

		run_pandoc ${kmer_alignment}
		mv $(basename ${kmer_alignment} | sed "s/.md/.html/g") ${output_dir}/${barcode_id}/de_novo_assembly_reads_subsets/reference_alignment/

	done

	cp ${barcode_dir}/ab1/*custom_consensus_assembly*.ab1 ${output_dir}/${barcode_id}/de_novo_assembly_reads_subsets/ab1/
	cp ${barcode_dir}/ab1/*trycycler_consensus_assembly*.ab1 ${output_dir}/${barcode_id}/de_novo_assembly_reads_subsets/ab1/

	# Run summary

	mv run_summaries_tmp/${barcode_id}_summary.tsv ${output_dir}/${barcode_id}/

	# Results files description

	cp docs/results_structure.html ${output_dir}/${barcode_id}/

done

rm -r run_summaries_tmp
