process FlyeAssembly {

  // Flye assembly
  
  label 'flye'

  publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.flye_output_dir}", mode: "copy", pattern: "*.fasta"
  publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.flye_output_dir}", mode: "copy", pattern: "*.graph.*"
  publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.flye_output_dir}", mode: "copy", pattern: "*.txt"

  input:
  tuple val(sample_id), val(genome_size), path(fastq)

  output:
  tuple val(sample_id), path("${sample_id}_*_flye.contigs.fasta"), emit: flye_assembly
  tuple val(sample_id), path("${sample_id}_*_flye.graph.gfa")
  tuple val(sample_id), path("${sample_id}_*_flye.graph.gv")
  tuple val(sample_id), path("${sample_id}_*_flye.info.txt"), emit: flye_report

  """
  output_prefix=${sample_id}_\$(basename ${fastq} | sed "s/.fq.gz//g" | sed "s/.fq//g" | sed "s/.fastq//g")_flye

  flye \
  --genome-size ${genome_size} \
  --threads \$SLURM_CPUS_ON_NODE \
  --out-dir . \
  ${params.flye_parameters} \
  ${params.flye_input_type} \
  ${fastq}

  mv assembly.fasta \${output_prefix}.contigs.fasta
  mv assembly_graph.gfa \${output_prefix}.graph.gfa
  mv assembly_graph.gv \${output_prefix}.graph.gv
  mv assembly_info.txt \${output_prefix}.info.txt
  """

}