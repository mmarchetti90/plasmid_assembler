/*
Trycycler assemblies consensus
*/

// ----------------Workflow---------------- //

include { TrycyclerCluster } from '../../modules/local/trycycler/trycycler_cluster.nf'
include { TrycyclerPreReconcile } from '../../modules/local/trycycler/trycycler_pre_reconcile.nf'
include { TrycyclerReconcile } from '../../modules/local/trycycler/trycycler_reconcile.nf'
include { TrycyclerMSA } from '../../modules/local/trycycler/trycycler_msa.nf'
include { TrycyclerPartition } from '../../modules/local/trycycler/trycycler_partition.nf'
include { TrycyclerConsensus } from '../../modules/local/trycycler/trycycler_consensus.nf'

workflow CONSENSUS_FROM_SUBSAMPLES_TRYCYCLER {

  take:
  concatenated_cleaned_fastq
  deconcatenated_assembly

  main:
  // CLUSTERING --------------------------- //

  // Join reads and assemblies channels
  concatenated_cleaned_fastq
    .join(deconcatenated_assembly.groupTuple(by: 0), by: 0, remainder: false)
    .set{ trycycler_cluster_input }

  // Cluster assemblies
  TrycyclerCluster(trycycler_cluster_input)

  trycycler_cluster = TrycyclerCluster.out.trycycler_cluster.transpose()

  // RECONCILING -------------------------- //

  // Pre-reconcile script channel
  pre_reconcile_script = Channel.fromPath("${projectDir}/scripts/pre_reconcile/pre_reconcile_filtering.py")

  // Filtering contigs before reconciling
  TrycyclerPreReconcile(pre_reconcile_script, trycycler_cluster)

  trycycler_cleaned_cluster = TrycyclerPreReconcile.out.trycycler_cleaned_cluster

  // Join reads and cluster channels
  concatenated_cleaned_fastq
    .join(trycycler_cleaned_cluster, by: 0, remainder: false)
    .set{ trycycler_reconcile_input }

  // Reconcile
  TrycyclerReconcile(trycycler_reconcile_input)

  trycycler_msa_input = TrycyclerReconcile.out.trycycler_cluster_reconciled
  
  // MULTIPLE SEQUENCE ALIGNMENT ---------- //

  TrycyclerMSA(trycycler_msa_input)

  trycycler_msa = TrycyclerMSA.out.trycycler_msa

  // READS PARTITION ---------------------- //

  // Join reads and cluster channels
  concatenated_cleaned_fastq
    .combine(trycycler_msa, by: 0)
    .set{ trycycler_partition_input }

  // Reads partition
  TrycyclerPartition(trycycler_partition_input)

  trycycler_partitioned = TrycyclerPartition.out.trycycler_partitioned

  // CONSENSUS ---------------------------- //

  TrycyclerConsensus(trycycler_partitioned)

  trycycler_consensus_data = TrycyclerConsensus.out.trycycler_consensus_data

  emit:
  trycycler_consensus_data

}