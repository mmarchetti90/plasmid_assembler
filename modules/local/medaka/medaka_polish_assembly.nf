process MedakaPolishAssembly {

  // Polish Flye or Canu assembly with Medaka
  
  label 'medaka'

  publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.consensus_output_dir}", mode: "copy", pattern: "*_allreads.polished.contigs.fasta"

  input:
  tuple val(sample_id), path(fastq), path(assembly)

  output:
  tuple val(sample_id), path("${sample_id}_allreads.polished.contigs.fasta"), emit: polished_assembly
  tuple val(sample_id), path("calls_to_draft.bam"), path("calls_to_draft.bam.bai"), emit: draft_alignment
  
  """
  medaka_consensus \
  -i ${fastq} \
  -d ${assembly} \
  -t \$SLURM_CPUS_ON_NODE \
  --bacteria \
  -o . \
  ${params.medaka_consensus_params}

  mv consensus.fasta ${sample_id}_allreads.polished.contigs.fasta
  """

}