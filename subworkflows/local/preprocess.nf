/*
Fastq concatenation, filtering, and subsampling
*/

// ----------------Workflow---------------- //

include { ConcatenateFastq } from '../../modules/local/concatenate_fastq/concatenate_fastq.nf'
include { Filtlong } from '../../modules/local/filtlong/filtlong.nf'
include { TrycyclerSubsample } from '../../modules/local/trycycler/trycycler_subsample.nf'

workflow PREPROCESS {

  take:
  fastq_input
  plasmid_size

  main:
  // CONCATENATE FASTQ -------------------- //

  ConcatenateFastq(fastq_input)

  concatenated_fastq = ConcatenateFastq.out.concatenated_fastq

  // CLEAN FASTQ -------------------------- //

  Filtlong(concatenated_fastq)

  concatenated_cleaned_fastq = Filtlong.out.concatenated_cleaned_fastq

  // SUBSAMPLE ---------------------------- //

  // Join plasmid_size and concatenated_cleaned_fastq channels
  plasmid_size
    .join(concatenated_cleaned_fastq, by: 0, remainder: false)
    .set{ subsampler_input }

  // Subsample
  TrycyclerSubsample(subsampler_input)

  // Transpose output to tuple of individual fastq files
  subsampled_fastqs = TrycyclerSubsample.out.subsampled_fastqs.transpose()

  emit:
  concatenated_cleaned_fastq
  subsampled_fastqs

}