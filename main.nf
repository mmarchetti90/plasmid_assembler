#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
Pipeline for plasmid assembly from ONT reads using Flye and Canu
*/

// ----------------Workflow---------------- //

include { PLASMID_ASSEMBLY } from './workflows/local/plasmid_assembly.nf'

workflow {

  // Run workflow
  PLASMID_ASSEMBLY()

}