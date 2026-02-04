---
geometry: margin=1.5cm
papersize: letter
sansfont:
fontsize: 12pt
urlcolor: blue
toc:
toc-depth: 4
---

# Output files details

## Overview of output files.

<pre>
<b>barcodeXX</b>
Main results directory.
│
├── <b>alignment_based_assembly</b>
│   This directory contains the reference-based assembly.
│   │
│   ├── <b>ab1</b>
│   │   │
│   │   └── <b>barcodeXX_reference_based_assembly.ab1</b>
│   │       ABIF file for the reference-based assembly.
│   │
│   ├── <b>barcodeXX_consensus.fasta</b>
│   │   Reference-based assembly sequence.
│   │
│   └── <b>minimap2_alignment</b>
│       │
│       └── <b>barcodeXX.alignment.log</b>
│           Stats of reads alignment to reference.
│
├── <b>barcodeXX_summary.tsv</b>
│   Overview of available assemblies.
│   (Can be imported in MS Excel)
│
├── <b>de_novo_assembly_all_reads</b>
│   This directory contains a de-novo assembly generated using all available ONT reads.
│   │
│   ├── <b>ab1</b>
│   │   │
│   │   └── <b>barcodeXX_all_reads_assembly.ab1</b>
│   │       ABIF file for the reference-based assembly.
│   │
│   ├── <b>assembly_annotations</b>
│   │   Plannotate assembly annotation.
│   │   │
│   │   ├── <b>barcodeXX_allreads_pLann.csv</b>
│   │   │   Overview of detected features.
│   │   │   (Can be imported in MS Excel)
│   │   │
│   │   ├── <b>barcodeXX_allreads_pLann.gbk</b>
│   │   │   GenBank annotation file.
│   │   │
│   │   └── <b>barcodeXX_allreads_pLann.html</b>
│   │       Graphical overview of the annotation.
│   │
│   ├── <b>barcodeXX_allreads.polished.contigs.fasta</b>
│   │   De-novo assembly sequence.
│   │
│   └── <b>reference_alignment</b>
│       │
│       ├── <b>barcodeXX_allreads_alignment.html</b>
│       │   Alignment of the assembly to the reference fasta using a custom alignment algorithm.
│       │
│       └── <b>barcodeXX_allreads_blastn.txt</b>
│           Alignment of the assembly to the reference fasta using BLASTN.
│
├── <b>de_novo_assembly_reads_subsets</b>
│   This directory contains de-novo assemblies generated using N subsets of the ONT reads (N = 12
│   as default).
│   Each reads subset is used to generate separate assemblies, which are then used to generate a
│   consensus using Trycycler or a custom algorithm.
│   │
│   ├── <b>ab1</b>
│   │   │
│   │   ├── <b>barcodeXX_custom_consensus_assembly.ab1</b>
│   │   │   ABIF file for the de-novo consensus assembly generated using a custom algorithm.
│   │   │
│   │   └── <b>barcodeXX_trycycler_consensus_assembly.ab1</b>
│   │       ABIF file for the de-novo consensus assembly generated using Trycycler.
│   │
│   ├── <b>assembly_annotations</b>
│   │   Plannotate assembly annotation.
│   │   │
│   │   ├── <b>barcodeXX_cluster_001_trycycler.consensus_pLann.csv</b>
│   │   │   Overview of detected features in the consensus assembly generated using Trycycler.
│   │   │   (Can be imported in MS Excel)
│   │   │
│   │   │
│   │   ├── <b>barcodeXX_cluster_001_trycycler.consensus_pLann.gbk</b>
│   │   │   GenBank annotation file for the consensus assembly generated using Trycycler.
│   │   │
│   │   ├── <b>barcodeXX_cluster_001_trycycler.consensus_pLann.html</b>
│   │   │   Graphical overview of the annotation for the consensus assembly generated using
│   │   │   Trycycler.
│   │   │
│   │   ├── <b>barcodeXX_custom.consensus_pLann.csv</b>
│   │   │   Overview of detected features in the consensus assembly generated using a custom
│   │   │   algorithm.
│   │   │   (Can be imported in MS Excel)
│   │   │
│   │   ├── <b>barcodeXX_custom.consensus_pLann.gbk</b>
│   │   │   GenBank annotation file for the consensus assembly generated using a custom algorithm.
│   │   │
│   │   └── <b>barcodeXX_custom.consensus_pLann.html</b>
│   │       Graphical overview of the annotation for the consensus assembly generated using a
│   │       custom algorithm.
│   │
│   ├── <b>final_assemblies</b>
│   │   This directory contains the final consensus assemblies.
│   │   │
│   │   ├── <b>barcodeXX_cluster_001_trycycler.consensus.polished.contigs.fasta</b>
│   │   │   Consensus assembly generated using Trycycler.
│   │   │
│   │   └── <b>barcodeXX_custom.consensus.polished.contigs.fasta</b>
│   │       Consensus assembly generated using a custom algorithm.
│   │
│   ├── <b>intermediate_assemblies</b>
│   │   │  This directory contains intermediate assemblies generated from ONT reads subsets.
│   │   │
│   │   ├── <b>barcodeXX_sample_01_flye.contigs.deconcatenated.fasta</b>
│   │   │
│   │   |   <b>...</b>
│   │   |   <b>...</b>
│   │   |   <b>...</b>
│   │   │
│   │   └── <b>barcodeXX_sample_N_flye.contigs.deconcatenated.fasta</b>
│   │
│   └── <b>reference_alignment</b>
│       │
│       ├── <b>barcodeXX_cluster_001_trycycler.consensus_alignment.html</b>
│       │   Alignment of the Trycycler consensus assembly to the reference fasta using a custom
│       │   alignment algorithm.
│       │
│       ├── <b>barcodeXX_cluster_001_trycycler.consensus_blastn.txt</b>
│       │   Alignment of the Trycycler consensus assembly to the reference fasta using BLASTN.
│       │
│       ├── <b>barcodeXX_custom.consensus_alignment.html</b>
│       │   Alignment of the custom algorithm consensus assembly to the reference fasta using a
│       │   custom alignment algorithm.
│       │
│       └── <b>barcodeXX_custom.consensus_blastn.txt</b>
│           Alignment of the custom algorithm consensus assembly to the reference fasta using
│           BLASTN.
│
└── <b>results_structure.html</b>
	This file.
</pre>

## Notes

* If a reference is provided, the reference-based assemblie is always available.

* De-novo assemblies can often fail, hence the multiple strategies.

* If files are missing, their automatic generation step likely failed.\
  Manual generation is possible on request.

* The generation of multiple de-novo assemblies from ONT reads subsets is usually the preferred\
  method. However, Trycycler is a tad too sensitive and could fail even when an optimal consensus\
  could be achieved. So, a custom algorithm was added to generate consensus assemblies\
  independently of Trycycler, by simply generating a multiple alignment of similar intermediate\
  assemblies. This custom algorithm will always generate a consensus assembly.

* The <b>barcodeXX_summary.tsv</b> file contains useful statistics of available assemblies, how\
  many contigs they contain, and their lenghts.
