/*
Assembly refinement, annotation, and comparison to reference
*/

// ----------------Workflow---------------- //

include { MedakaPolishAssembly } from '../../modules/local/medaka/medaka_polish_assembly.nf'
include { SamtoolsPileup } from '../../modules/local/samtools/samtools_pileup.nf'
include { PileupToAb1 } from '../../modules/local/ab1/pileup_to_ab1.nf'
include { AnnotateAssembly } from '../../modules/local/plannotate/plannotate.nf'
include { PairedBlastNAlignment } from '../../modules/local/alignment/paired_blastn_alignment.nf'
include { KmerAlignment } from '../../modules/local/alignment/paired_kmer_alignment.nf'

workflow DOWNSTREAM_PROCESSING_FULL_SAMPLES {

  take:
  deconcatenated_assembly_full_samples
  concatenated_cleaned_fastq
  reference_fasta

  main:
  // POLISHING ASSEMBLY ------------------- //

  // Join fastq and assembly channels
  concatenated_cleaned_fastq
    .join(deconcatenated_assembly_full_samples, by: 0, remainder: false)
    .set{ medaka_input }

  // Polish assembly
  MedakaPolishAssembly(medaka_input)

  polished_assembly = MedakaPolishAssembly.out.polished_assembly
  draft_alignment = MedakaPolishAssembly.out.draft_alignment

  // AB1-LIKE FILE GENERATION ------------- //

  // Samtools pileup for draft_alignment
  SamtoolsPileup(draft_alignment)

  alignment_pileup = SamtoolsPileup.out.alignment_pileup

  // ab1 generation script channel
  alignment_script = Channel.fromPath("${projectDir}/scripts/ab1/ab1writer.py")

  // ab1 generation
  PileupToAb1(alignment_script, "all_reads_assembly", alignment_pileup)

  // ANNOTATION --------------------------- //

  AnnotateAssembly(polished_assembly)

  // BLASTN ALIGNMENT --------------------- //

  // Join reference and assembly channels
  reference_fasta
    .join(polished_assembly, by: 0, remainder: false)
    .set{ reference_and_assembly }

  // BlastN alignment
  PairedBlastNAlignment(reference_and_assembly)

  // KMER ALIGNMENT ----------------------- //

  // Kmer alignment script channel
  alignment_script = Channel.fromPath("${projectDir}/scripts/alignment/kmer_alignment.py")

  // Kmer alignment
  KmerAlignment(alignment_script, reference_and_assembly)

}