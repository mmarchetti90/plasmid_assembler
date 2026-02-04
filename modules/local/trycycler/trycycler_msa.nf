process TrycyclerMSA {

  // Assemblies MSA with Trycycler
  
  label 'trycycler'

  input:
  tuple val(sample_id), path(cluster_dir)

  output:
  tuple val(sample_id), path("cluster_*_msa"), emit: trycycler_msa
  
  """
  cluster_id=\$(basename ${cluster_dir} | sed "s/_reconciled//g")

  cp -r -L ${cluster_dir} \${cluster_id}_msa

  trycycler msa \
  --threads \$SLURM_CPUS_ON_NODE \
  --cluster_dir \${cluster_id}_msa \
  ${params.trycycler_msa_params}
  """

}