/*
Plasmid assembly from ONT reads using Flye or Canu and subsampled reads
*/

// ----------------Workflow---------------- //

include { CanuAssembly } from '../../modules/local/canu/canu_assembly.nf'
include { FlyeAssembly } from '../../modules/local/flye/flye_assembly.nf'
include { DeconcatenateAssemblies } from '../../modules/local/deconcatenate_assemblies/deconcatenate_assemblies.nf'

workflow ASSEMBLE_FULL_SAMPLES {

  take:
  concatenated_cleaned_fastq
  plasmid_size

  main:
  // ASSEMBLY ----------------------------- //

  // Combine plasmid_size and concatenated_cleaned_fastq channels
  plasmid_size
    .join(concatenated_cleaned_fastq, by: 0, remainder: false)
    .set{ assembler_input }

  if (params.assembler == "flye") {

    // Flye assembly
    FlyeAssembly(assembler_input)

    assembly = FlyeAssembly.out.flye_assembly

  }
  else {

    // Canu assembly
    CanuAssembly(assembler_input)

    assembly = CanuAssembly.out.canu_assembly

  }

  // DECONCATENATE ASSEMBLIES ------------- //

  // Deconcatenation script channel
  deconcatenation_script = Channel.fromPath("${projectDir}/scripts/deconcatenation/deconcatenate_assemblies.py")

  // Join plasmid_size and assembly channels
  plasmid_size
    .join(assembly, by: 0, remainder: false)
    .set{ plasmid_size_and_assembly }

  // Assembly deconcatenation
  DeconcatenateAssemblies(deconcatenation_script, plasmid_size_and_assembly)

  deconcatenated_assembly = DeconcatenateAssemblies.out.deconcatenated_assembly

  emit:
  deconcatenated_assembly

}