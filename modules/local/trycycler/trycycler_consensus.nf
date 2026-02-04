process TrycyclerConsensus {

  // Assemblies consensus with Trycycler
  
  label 'trycycler'

  publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.trycycler_output_dir}", mode: "copy", pattern: "*_cluster_*_trycycler"

  input:
  tuple val(sample_id), path(cluster_dir)

  output:
  tuple val(sample_id), path("${sample_id}_cluster_*_trycycler"), emit: trycycler_consensus_data
  
  """
  cluster_id=\$(basename ${cluster_dir} | sed "s/_partitioned//g")

  cp -r -L ${cluster_dir} ${sample_id}_\${cluster_id}_trycycler

  trycycler consensus \
  --threads \$SLURM_CPUS_ON_NODE \
  --cluster_dir ${sample_id}_\${cluster_id}_trycycler \
  ${params.trycycler_consensus_params}
  """

}