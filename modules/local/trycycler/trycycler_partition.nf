process TrycyclerPartition {

  // Reads partition with Trycycler
  
  label 'trycycler'

  //publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.trycycler_output_dir}", mode: "copy", pattern: "*_assemblies_cluster_reconciled_msa_partitioned"

  input:
  tuple val(sample_id), path(fastq), path(cluster_dir)

  output:
  tuple val(sample_id), path("cluster_*_partitioned"), emit: trycycler_partitioned

  """
  cluster_id=\$(basename ${cluster_dir} | sed "s/_msa//g")

  cp -r -L ${cluster_dir} \${cluster_id}_partitioned

  trycycler partition \
  --threads \$SLURM_CPUS_ON_NODE \
  --reads ${fastq} \
  --cluster_dirs \${cluster_id}_partitioned \
  ${params.trycycler_partition_params}
  """

}