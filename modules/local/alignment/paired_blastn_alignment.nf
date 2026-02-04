process PairedBlastNAlignment {

  // Align assembly to reference
  
  label 'blast'

  publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.blastn_output_dir}", mode: "copy", pattern: "*_blastn.txt"

  input:
  tuple val(sample_id), path(reference), path(assembly)

  output:
  tuple val(sample_id), path("*_blastn.txt"), emit: blast_alignment
  
  """
  sample_prefix=\$(basename ${assembly} | \
  sed "s/.polished//g" | \
  sed "s/.contigs.fasta//g")

  blastn \
  -query ${reference} \
  -subject ${assembly} \
  -outfmt 0 \
  -out \${sample_prefix}_blastn.txt
  """

}