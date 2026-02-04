process TrycyclerCluster {

  // Assemblies clustering with Trycycler
  
  label 'trycycler'

  //publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.trycycler_output_dir}", mode: "copy", pattern: "*_assemblies_cluster"

  input:
  tuple val(sample_id), path(fastq), path(assembly)

  output:
  tuple val(sample_id), path("trycycler_clusters/cluster_*"), emit: trycycler_cluster
  tuple val(sample_id), path("${sample_id}_contigs.phylip")
  tuple val(sample_id), path("${sample_id}_contigs.newick")

  """
  trycycler cluster \
  --threads \$SLURM_CPUS_ON_NODE \
  --assemblies *.fasta \
  --reads ${fastq} \
  --out_dir trycycler_clusters \
  ${params.trycycler_cluster_params}

  mv trycycler_clusters/contigs.phylip ${sample_id}_contigs.phylip
  mv trycycler_clusters/contigs.newick ${sample_id}_contigs.newick
  """

}