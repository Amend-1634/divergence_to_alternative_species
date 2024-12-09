# Usage

This Nextflow pipeline is designed for remapping the reads to other fasta files and estimate the number of mismatches, after removing all potential damage (C-to-T at forward strand, G-to-A at reverse strand).


## Setup

This pipeline has been tested with Nextflow version 22.10.1. Please ensure you have the appropriate version installed before running the pipeline.

## To run the workflow


nextflow run tip_dating_1.nf \
        --label "$(basename $file | sed 's/.fas//')" \
        --all_input "/path/to/genus.fas" \
        --threads "5" \
        -resume

## Input

example for ${taxa}.fas: first row as the path to the bam file with all mapped reads, the second row as the reference for the bam file, and the rest of the row are other alternative reference genomes for the mapped reads.

```
$ cat Poa.fas 
/path/to/61genus-Poa_pratensis_4545-ext.bam
/path/to/Poa_pratensis_4545.fna
/path/to/Poa_pratensis_subsp._pratensis_368382.fna
/path/to/Poa_glauca_227214.fna
```

## Usage

This is based on the assumption that closely related taxa will have reads mapped to both of them and we aim at inferring which is the most likely species by account for the number of mutations divided by the reference genome size. This approach supported by individual mapping is to not lose mutation information for ambiguiously mapped reads in competitive mapping.
