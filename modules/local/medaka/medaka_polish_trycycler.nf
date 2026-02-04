process MedakaPolishTrycycler {

  // Polish trycycler consensus assembly with Medaka
  
  label 'medaka'

  publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.consensus_output_dir}", mode: "copy", pattern: "*_trycycler.consensus.polished.contigs.fasta"

  input:
  tuple val(sample_id), path(cluster_dir)

  output:
  tuple val(sample_id), path("${sample_id}_cluster_*_trycycler.consensus.polished.contigs.fasta"), emit: polished_assembly
  tuple val(sample_id), path("calls_to_draft.bam"), path("calls_to_draft.bam.bai"), emit: draft_alignment
  
  """
  cluster_id=\$(basename ${cluster_dir} | sed "s/_trycycler//g" | sed "s/${sample_id}_//g")

  medaka_consensus \
  -i ${cluster_dir}/4_reads.fastq \
  -d ${cluster_dir}/7_final_consensus.fasta \
  -t \$SLURM_CPUS_ON_NODE \
  --bacteria \
  -o . \
  ${params.medaka_consensus_params}

  mv consensus.fasta ${sample_id}_\${cluster_id}_trycycler.consensus.polished.contigs.fasta
  """

}