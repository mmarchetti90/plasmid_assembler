process TrycyclerSubsample {

  // Subset fastq
  
  label 'trycycler'

  //publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.concatenated_fastq_dir}", mode: "copy", pattern: "*.fq"

  input:
  tuple val(sample_id), val(genome_size), path(fastq)

  output:
  tuple val(sample_id), path("read_subsets/*.{fq,fastq}"), emit: subsampled_fastqs

  """
  trycycler subsample \
  --threads \$SLURM_CPUS_ON_NODE \
  --genome_size ${genome_size} \
  ${params.subsampling_parameters} \
  --reads ${fastq} \
  --out_dir read_subsets
  """

}