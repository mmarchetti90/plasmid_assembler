process SamtoolsConsensus {

  // Runs samtools consensus
  
  label 'samtools'

  publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.alignment_based_assembly_dir}", mode: "copy", pattern: "*_consensus.fasta"

  input:
  tuple val(sample_id), path(bam), path(bai)

  output:
  path "*_consensus.fasta", optional: false, emit: alignment_consensus
  
  """
  samtools consensus \
  -@ \$SLURM_CPUS_ON_NODE \
  --show-del yes \
  --show-ins yes \
  --format fasta \
  ${bam} > ${sample_id}_consensus.fasta
  """

}