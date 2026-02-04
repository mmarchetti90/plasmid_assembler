process SamtoolsPileup {

  // Runs samtools pileup
  
  label 'samtools'

  //publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.alignment_based_assembly_dir}", mode: "copy", pattern: "*_pileup.txt"

  input:
  tuple val(sample_id), path(bam), path(bai)

  output:
  tuple val(sample_id), path("*_pileup.txt"), optional: false, emit: alignment_pileup
  
  """
  samtools mpileup \
  ${bam} > ${sample_id}_pileup.txt
  """

}