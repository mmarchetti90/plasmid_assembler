/*
Workflow for plasmid assembly from ONT reads using Flye and Canu
*/

// ----------------Workflow---------------- //

include { PREPROCESS } from '../../subworkflows/local/preprocess.nf'
include { ASSEMBLE_FULL_SAMPLES } from '../../subworkflows/local/assemble_full_samples.nf'
include { ASSEMBLE_SUBSAMPLES } from '../../subworkflows/local/assemble_subsamples.nf'
include { CONSENSUS_FROM_SUBSAMPLES_TRYCYCLER } from '../../subworkflows/local/consensus_from_subsamples_trycycler.nf'
include { CONSENSUS_FROM_SUBSAMPLES_CUSTOM } from '../../subworkflows/local/consensus_from_subsamples_custom.nf'
include { DOWNSTREAM_PROCESSING_FULL_SAMPLES } from '../../subworkflows/local/downstream_processing_full_samples.nf'
include { DOWNSTREAM_PROCESSING_SUBSAMPLES_TRYCYCLER } from '../../subworkflows/local/downstream_processing_subsamples_trycycler.nf'
include { DOWNSTREAM_PROCESSING_SUBSAMPLES_CUSTOM } from '../../subworkflows/local/downstream_processing_subsamples_custom.nf'
include { REFERENCE_BASED_ASSEMBLY } from '../../subworkflows/local/reference_based_assembly.nf'

workflow PLASMID_ASSEMBLY {

  main:
  // INPUT CHANNELS ----------------------- //

  Channel
    .fromPath("${params.sample_manifest_path}")
    .splitCsv(header: true, sep: '\t')
    .map{ row -> tuple(row.sample_id, file(row.sample_path)) }
    .set{ fastq_input }

  Channel
    .fromPath("${params.sample_manifest_path}")
    .splitCsv(header: true, sep: '\t')
    .map{ row -> tuple(row.sample_id, row.plasmid_size) }
    .set{ plasmid_size }

  Channel
    .fromPath("${params.sample_manifest_path}")
    .splitCsv(header: true, sep: '\t')
    .map{ row -> tuple(row.sample_id, file(row.reference_fasta)) }
    .set{ reference_fasta }

  // PREPROCESSING ------------------------ //

  PREPROCESS(fastq_input, plasmid_size)

  concatenated_cleaned_fastq = PREPROCESS.out.concatenated_cleaned_fastq
  subsampled_fastqs = PREPROCESS.out.subsampled_fastqs

  // ASSEMBLY ----------------------------- //

  // Assembly of full samples
  ASSEMBLE_FULL_SAMPLES(concatenated_cleaned_fastq, plasmid_size)

  deconcatenated_assembly_full_samples = ASSEMBLE_FULL_SAMPLES.out.deconcatenated_assembly

  // Assembly of subsamples
  ASSEMBLE_SUBSAMPLES(subsampled_fastqs, plasmid_size)

  deconcatenated_assembly_subsamples = ASSEMBLE_SUBSAMPLES.out.deconcatenated_assembly

  // CONSENSUS FROM SUBSAMPLES ------------ //

  // Trycycler
  CONSENSUS_FROM_SUBSAMPLES_TRYCYCLER(concatenated_cleaned_fastq, deconcatenated_assembly_subsamples)

  trycycler_consensus_data = CONSENSUS_FROM_SUBSAMPLES_TRYCYCLER.out.trycycler_consensus_data

  // Custom implementation
  CONSENSUS_FROM_SUBSAMPLES_CUSTOM(deconcatenated_assembly_subsamples)

  custom_consensus_fasta = CONSENSUS_FROM_SUBSAMPLES_CUSTOM.out.consensus_fasta

  // DOWNSTREAM PROCESSING ---------------- //

  // Full samples assemblies
  DOWNSTREAM_PROCESSING_FULL_SAMPLES(deconcatenated_assembly_full_samples, concatenated_cleaned_fastq, reference_fasta)

  // Trycycler consensus assemblies
  DOWNSTREAM_PROCESSING_SUBSAMPLES_TRYCYCLER(trycycler_consensus_data, reference_fasta)

  // Custom consensus assemblies
  DOWNSTREAM_PROCESSING_SUBSAMPLES_CUSTOM(custom_consensus_fasta, concatenated_cleaned_fastq, reference_fasta)

  // REFERENCE-BASED ASSEMBLY ------------- //

  REFERENCE_BASED_ASSEMBLY(reference_fasta, concatenated_cleaned_fastq)

}