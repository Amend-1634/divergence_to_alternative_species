# Tip Dating Pipeline with 3bp End Ignoring

This Nextflow pipeline is designed for estimating its sequence divergence to all alternative reference taxa, specifically ignoring the last 3 base pairs (3bp) on either side of the reads. This is to avoid bias in the estimation caused by ancient DNA damage at the read ends, following the particular protocol applied to the Alazeya steppe bison gut content sample.

## Setup

This pipeline has been tested with Nextflow version 22.10.1. Please ensure you have the appropriate version installed before running the pipeline.

## To run the workflow

while IFS= read -r genus;do
        echo ${dir}${genus}.fas

        time nextflow run tip_dating_3bp.nf \
        --label "$genus" \
        --all_input "fas/${genus}.fas" \
        --threads "15" \
        -resume

done < list_genus

## Input

example for ${taxa}.fas: first row as the path to the bam file with all mapped reads, the second row as the reference for the bam file, and the rest of the row are other alternative reference genomes for the mapped reads.

```
$ cat Poa.fas 
/crex/proj/snic2022-6-144/nobackup/CHENYU/tip_dating/data1/bam/61genus-Poa_pratensis_4545-ext.bam
/crex/proj/snic2022-6-144/nobackup/CHENYU/tip_dating/data1/fasta/Poa_pratensis_4545.fna
/crex/proj/snic2022-6-144/nobackup/CHENYU/tip_dating/data1/fasta/Poa_pratensis_subsp._pratensis_368382.fna
/crex/proj/snic2022-6-144/nobackup/CHENYU/tip_dating/data1/fasta/Poa_glauca_227214.fna
```

## Usage

This is based on the assumption that closely related taxa will have reads mapped to both of them and we aim at inferring which is the most likely species by account for the number of mutations divided by the reference genome size. This approach supported by individual mapping is to not lose mutation information for ambiguiously mapped reads in competitive mapping.
