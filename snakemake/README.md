### Description

This folder contains snakemake variant calling workflow.


### Workflow details

The main workflow is in `variants.smk` file.
This workflow includes other workflows such as `bwa.smk` and `bcftools.smk`.

`bwa.smk` : Alignment workflow for both single-end and paired-end reads using bwa-mem program.

`bcftools.smk` : Variant calling workflow using bcftools.

`variant.smk` : Master variant calling workflow that combines bwa.smk and bcftools.smk

### How to run?

Run using 4 processors and parameters specified in config file.

    snakemake -s variants.smk -c 4

To override any of the runtime parameters specified in the config file, for eg: alignment directory

    snakemake -s variants.smk -c 4 --config aln_dir=align

To perform a dry-run

    snakemake -s variants.smk -c 4 --config aln_dir=align -n

### Sample sheet specification

Sample sheet is a comma-separated file with the header sample,read1,read2

An example sample sheet would look like this:

    sample,read1,read2
    S1,data/S1_R1.fq,data/S1_R2.fq
    S2,data/S2_R1.fq,data/S2_R2.fq


### Workflow requirements

1. Software requirements 
   1. Snakemake
   2. bwa
   3. bcftools
2. Configuration requirements
   1. All parameters are specified through a YAML configuration file, eg: config.yaml 
   2. A samplesheet in comma-separated format with atleast 2 columns with column headers as 'sample,read1,read2'


An example sample sheet would be

    sample,read1,read2
    S1,data/S1_R1.fq,data/S1_R2.fq
    S2,data/S2_R1.fq,data/S2_R2.fq`

