process KmerAlignment {

  // Align assembly to reference
  
  label 'python'

  publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.kmer_alignment_output_dir}", mode: "copy", pattern: "*_alignment.md"

  input:
  each path(alignment_script)
  tuple val(sample_id), path(reference), path(assembly)

  output:
  tuple val(sample_id), path("*_alignment.md"), emit: blast_alignment
  
  """
  sample_prefix=\$(basename ${assembly} | \
  sed "s/.polished//g" | \
  sed "s/.contigs.fasta//g")

  python ${alignment_script} \
  --query ${reference} \
  --subject ${assembly}

  mv alignment.md \${sample_prefix}_alignment.md
  """

}