process ConcatenateFastq {

  // Concatenate fastq files into one
  
  label 'local'

  //publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.concatenated_fastq_dir}", mode: "copy", pattern: "*.fq.gz"

  input:
  tuple val(sample_id), path(sample_path)

  output:
  tuple val(sample_id), path("${sample_id}_allreads.fq.gz"), emit: concatenated_fastq

  """
  if [[ "${sample_path}" == *.gz ]]
  then

    cp ${sample_path} ${sample_id}_allreads.fq.gz

  else

    zcat ${sample_path}/*.gz | gzip > ${sample_id}_allreads.fq.gz

  fi
  """

}