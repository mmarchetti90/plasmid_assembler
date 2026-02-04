process PileupToAb1 {

  // Converts a samtools mpileup file to ab1-like trace
  
  label 'python'

  publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.ab1_output_dir}", mode: "copy", pattern: "*.ab1"

  input:
  each path(abiwriter_script)
  val pileup_type
  tuple val(sample_id), path(pileup)

  output:
  tuple val(sample_id), path("${sample_id}_${pileup_type}.ab1"), emit: synthetic_ab1_file
  
  """
  python ${abiwriter_script} \
  --pileup ${pileup} \
  --output ${sample_id}_${pileup_type}.ab1
  """

}