process SamtoolsFlagstat {

  // Runs samtools flagstat
  
  label 'samtools'

  publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.minimap2_alignment_dir}", mode: "copy", pattern: "*.alignment.log"

  input:
  tuple val(sample_id), path(bam), path(bai)

  output:
  path "*.alignment.log", optional: false, emit: alignment_reports
  
  """
  samtools flagstat \
  -@ \$SLURM_CPUS_ON_NODE \
  -O 'default' \
  ${bam} > ${sample_id}.alignment.log
  """

}