/*
Workflow for plasmid assembly from ONT reads using a reference fasta
*/

// ----------------Workflow---------------- //

include { Minimap2Alignment } from '../../modules/local/alignment/minimap2_alignment.nf'
include { SamToBam } from '../../modules/local/samtools/sam_to_bam.nf'
include { SamtoolsFlagstat } from '../../modules/local/samtools/samtools_flagstat.nf'
include { SamtoolsConsensus } from '../../modules/local/samtools/samtools_consensus.nf'
include { SamtoolsPileup } from '../../modules/local/samtools/samtools_pileup.nf'
include { PileupToAb1 } from '../../modules/local/ab1/pileup_to_ab1.nf'

workflow REFERENCE_BASED_ASSEMBLY {

  take:
  reference_fasta
  concatenated_cleaned_fastq

  main:
  // MINIMAP2 ON REFERENCE ---------------- //

  // Merge reference_fasta and concatenated_cleaned_fastq channels
  reference_fasta
    .join(concatenated_cleaned_fastq, by: 0, remainder: false)
    .set{ minimap2_reference_input }

  // Minimap2 on reference fasta
  Minimap2Alignment(minimap2_reference_input)

  reference_aligned_sam = Minimap2Alignment.out.aligned_sam

  // Sam to bam from reference alignment
  SamToBam(reference_aligned_sam)

  reference_aligned_sorted_bam = SamToBam.out.aligned_sorted_bam

  // Reference alignment stats
  SamtoolsFlagstat(reference_aligned_sorted_bam)

  // Samtools consensus for reference aligned bam
  SamtoolsConsensus(reference_aligned_sorted_bam)

  // AB1-LIKE FILE GENERATION ------------- //

  // Samtools pileup for reference aligned bam
  SamtoolsPileup(reference_aligned_sorted_bam)

  reference_alignment_pileup = SamtoolsPileup.out.alignment_pileup

  // ab1 generation script channel
  alignment_script = Channel.fromPath("${projectDir}/scripts/ab1/ab1writer.py")

  // ab1 generation
  PileupToAb1(alignment_script, "reference_based_assembly", reference_alignment_pileup)

}