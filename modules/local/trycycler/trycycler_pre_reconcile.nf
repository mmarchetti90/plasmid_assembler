process TrycyclerPreReconcile {

  // Cleaning Trycycler cluster before reconciling
  // Reduces the chances of trycycler reconcile failing
  
  label 'python'

  input:
  each path(pre_reconcile_script)
  tuple val(sample_id), path(cluster_dir)

  output:
  tuple val(sample_id), path("cluster_*_cleaned"), emit: trycycler_cleaned_cluster

  """
  cluster_id=\$(basename ${cluster_dir})

  assemblies=tmp
  
  for fasta in ${cluster_dir}/1_contigs/*.fasta
  do
  
    assemblies=\${assemblies},\${fasta}
  
  done
  
  assemblies=\$(echo \${assemblies} | sed "s/tmp,//g")

  python ${pre_reconcile_script} \
  --assemblies \${assemblies}

  mkdir -p \${cluster_id}_cleaned/1_contigs
  mv *.fasta \${cluster_id}_cleaned/1_contigs/
  """

}