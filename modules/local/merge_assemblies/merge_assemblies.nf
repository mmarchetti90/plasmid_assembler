process CustomConsensus {

  // Assemblies consensus using a custom script
  
  label 'python'

  publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.consensus_output_dir}", mode: "copy", pattern: "*_custom.consensus.contigs.fasta"
  publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.consensus_output_dir}", mode: "copy", pattern: "*_custom.subassemblies.alignments.txt"

  input:
  each path(merging_script)
  tuple val(sample_id), path(assembly)

  output:
  tuple val(sample_id), path("*_custom.consensus.contigs.fasta"), emit: consensus_fasta
  tuple val(sample_id), path("*_custom.subassemblies.alignments.txt"), emit: subassemblies_alignment

  """
  python ${merging_script} \
  --sample_id ${sample_id} \
  --fasta_dir .
  """

}