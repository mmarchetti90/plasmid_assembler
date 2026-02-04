process CanuAssembly {

  // Canu assembly
  
  label 'canu'

  publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.canu_output_dir}", mode: "copy", pattern: "*.report"
  publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.canu_output_dir}", mode: "copy", pattern: "*.fastq.gz"
  publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.canu_output_dir}", mode: "copy", pattern: "*.fasta"
  publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.canu_output_dir}", mode: "copy", pattern: "*.layout.*"

  input:
  tuple val(sample_id), val(genome_size), path(fastq)

  output:
  tuple val(sample_id), path("${sample_id}_*_canu.report"), emit: canu_report
  tuple val(sample_id), path("${sample_id}_*_canu.correctedReads.fasta.gz")
  tuple val(sample_id), path("${sample_id}_*_canu.trimmedReads.fasta.gz")
  tuple val(sample_id), path("${sample_id}_*_canu.contigs.fasta"), emit: canu_assembly
  tuple val(sample_id), path("${sample_id}_*_canu.unassembled.fasta")
  tuple val(sample_id), path("${sample_id}_*_canu.contigs.layout.readToTig")
  tuple val(sample_id), path("${sample_id}_*_canu.contigs.layout.tigInfo")

  """
  output_prefix=${sample_id}_\$(basename ${fastq} | sed "s/.fq.gz//g" | sed "s/.fq//g" | sed "s/.fastq//g")_canu

  canu \
  -p \${output_prefix} \
  -maxThreads=\$SLURM_CPUS_ON_NODE \
  genomeSize=${genome_size} \
  -nanopore \
  ${params.canu_parameters} \
  ${fastq}
  """

}