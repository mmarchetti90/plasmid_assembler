process TrycyclerReconcile {

  // Assemblies reconciling with Trycycler
  
  label 'trycycler'

  input:
  tuple val(sample_id), path(fastq), path(cluster_dir)

  output:
  tuple val(sample_id), path("cluster_*_reconciled"), emit: trycycler_cluster_reconciled

  """
  cluster_id=\$(basename ${cluster_dir} | sed "s/_cleaned//g")

  cp -r -L ${cluster_dir} \${cluster_id}_reconciled

  trycycler reconcile \
  --threads \$SLURM_CPUS_ON_NODE \
  --reads ${fastq} \
  --cluster_dir \${cluster_id}_reconciled \
  ${params.trycycler_reconcile_params}
  """

}