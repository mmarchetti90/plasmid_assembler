/*
Assemblies consensus using a custom script
*/

// ----------------Workflow---------------- //

include { CustomConsensus } from '../../modules/local/merge_assemblies/merge_assemblies.nf'

workflow CONSENSUS_FROM_SUBSAMPLES_CUSTOM {

  take:
  deconcatenated_assembly

  main:
  // CONSENSUS ---------------------------- //

  // Assemblies merging script channel
  merging_script = Channel.fromPath("${projectDir}/scripts/merging/merge_assemblies.py")

  // Group assemblies by sample id
  deconcatenated_assembly.groupTuple(by: 0)
    .set{ grouped_assemblies }

  // Consensus
  CustomConsensus(merging_script, grouped_assemblies)

  consensus_fasta = CustomConsensus.out.consensus_fasta

  emit:
  consensus_fasta

}