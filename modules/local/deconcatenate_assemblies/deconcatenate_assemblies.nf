process DeconcatenateAssemblies {

  // Deconcatenate assemblies
  
  label 'python'

  publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.deconcatenated_output_dir}", mode: "copy", pattern: "*.log"
  publishDir "${projectDir}/${params.main_output_dir}/${sample_id}/${params.deconcatenated_output_dir}", mode: "copy", pattern: "*.fasta"

  input:
  each path(deconcatenation_script)
  tuple val(sample_id), val(genome_size), path(assembly)

  output:
  tuple val(sample_id), path("*.contigs.deconcatenated.fasta"), emit: deconcatenated_assembly

  """
  output_prefix=\$(basename ${assembly} | sed "s/.contigs.fasta//g")

  python ${deconcatenation_script} \
  --assembly ${assembly} \
  --expected_size ${genome_size} &> \${output_prefix}_deconcatenation.log

  mv deconcatenated.fasta \${output_prefix}.contigs.deconcatenated.fasta
  """

}