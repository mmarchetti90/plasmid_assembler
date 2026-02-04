process Filtlong {

  // Filter fastq
  
  label 'filtlong'

  //publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.concatenated_fastq_dir}", mode: "copy", pattern: "*.fq"

  input:
  tuple val(sample_id), path(fastq)

  output:
  tuple val(sample_id), path("${sample_id}_allreads.cleaned.fq"), emit: concatenated_cleaned_fastq

  """
  filtlong \
  ${params.filtlong_parameters} \
  ${fastq} > ${sample_id}_allreads.cleaned.fq
  """

}