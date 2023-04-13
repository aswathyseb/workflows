### Description

This folder contains snakemake variant calling workflow.


### Workflow details

The main workflow is in `variants.smk` file.
This workflow includes other workflows such as `bwa.smk` and `bcftools.smk`.

bwa.smk : Alignment workflow for both single-end and paired-end reads using bwa-mem program
bcftools.smk : Variant calling workflow using bcftools.
variant.smk : Master variant calling workflow that combines bwa.smk and bcftools.smk


### Workflow requirements

1. Software requirements 
   1. Snakemake
   2. bwa
   3. bcftools
2. Configuration requirements
   1. All parameters are specified through a YAML configuration file, eg: config.yaml 
   2. A samplesheet in tab-delimited format with atleast 2 columns with column headers as 'sample\tread1\tread2'


An example sample sheet would be

    sample	read1	read2
    S1	data/S1_R1.fq	data/S1_R2.fq
    S2	data/S2_R1.fq	data/S2_R2.fq`

