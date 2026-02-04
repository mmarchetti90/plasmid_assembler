process AnnotateAssembly {

  // Annotate assembly using pLannotate
  
  label 'plannotate'

  publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.annotated_assembly_dir}", mode: "copy", pattern: "*_pLann.{gbk,csv,html}"

  input:
  tuple val(sample_id), path(assembly)

  output:
  tuple val(sample_id), path("*_pLann.gbk"), emit: plannotate_gbk
  tuple val(sample_id), path("*_pLann.csv"), emit: plannotate_csv
  tuple val(sample_id), path("*_pLann.html"), emit: plannotate_html
  
  """
  sample_prefix=\$(basename ${assembly} | \
  sed "s/.polished//g" | \
  sed "s/.contigs.fasta//g")

  plannotate batch \
  -i ${assembly} \
  -o . \
  -f \${sample_prefix} \
  --csv \
  --html
  """

}