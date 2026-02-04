process Minimap2Alignment {

  // Align fastq to reference
  
  label 'minimap2'

  //publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.minimap2_alignment_dir}", mode: "copy", pattern: "*.aligned.sam"

  input:
  tuple val(sample_id), path(reference), path(ont_reads)

  output:
  tuple val("${sample_id}"), path("${sample_id}.aligned.sam"), optional: false, emit: aligned_sam
  
  """
  minimap2 \
  ${params.minimap2_params} \
  -t \$SLURM_CPUS_ON_NODE \
  -ax map-ont \
  ${reference} \
  ${ont_reads} \
  -y \
  -a > ${sample_id}.aligned.sam
  """

}