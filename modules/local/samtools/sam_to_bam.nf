process SamToBam {

  // Converts a sam file to a sorted, indexed bam
  
  label 'samtools'

  publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.minimap2_alignment_dir}", mode: "copy", pattern: "*.aligned.sorted.{bam,bam.bai}"

  input:
  tuple val(sample_id), path(sam)

  output:
  tuple val("${sample_id}"), path("${sample_id}.aligned.sorted.bam"), path("${sample_id}.aligned.sorted.bam.bai"), optional: false, emit: aligned_sorted_bam
  
  """
  # Sort and convert to bam
  samtools sort \
  -@ \$SLURM_CPUS_ON_NODE \
  ${sam} \
  -o ${sample_id}.aligned.sorted.bam

  # Index
  samtools index \
  -b \
  -@ \$SLURM_CPUS_ON_NODE \
  ${sample_id}.aligned.sorted.bam
  """

}