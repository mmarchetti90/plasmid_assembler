# De-novo plasmid assembly from ONT reads

## Overview

The pipeline attempts to assemble sequences using 3 approaches: using all samples's sequences; using subsampled sequences; using a fasta reference.

### Preprocessing

* Concatenation of fastq files from the same sample.
* Fastq quality filtering.
* Fastq subsampling.

### Assembly using all sequences

* Assembly with Flye or Canu.
* De-concatenation.
* Assembly refinement with Medaka.
* Comparison to reference sequence (if fasta is provided).
* Annotation with pLannotate.
* ABIF trace generation.

### Assembly using fastq subsamples

* Assembly with Flye or Canu.
* De-concatenation.
* Consensus assembly with Trycycler and a custom algorithm, in parallel.
* Assembly refinement with Medaka.
* Comparison to reference sequence (if fasta is provided).
* Annotation with pLannotate.
* ABIF trace generation.

### Reference-based assembly

**N.B.** Optional, only run if a reference is provided.

* Minimap2 alignment.
* Consensus fasta generation from bam.
* ABIF trace generation.

## Input

The pipeline reads in a tab-separated file with the following required columns:

* **sample_id** : unique sample identifier that will be used for naming outputs.
* **sample_path** : path to either a fastq file or a folder with fastq files to be concatenated.
* **genome_size** : approximate size of the plasmid in base-pairs.
* **reference_fasta** : path to reference fasta for the plasmid. If missing, set to "mock.fasta".

Note that input fastq files **must** be be gzipped.

## Output

See docs/results_structure.{md,html} for an in-depth description.

## Notes

* Flye and Canu were used because they are popular, well-maintained, and easy to install (conda).
