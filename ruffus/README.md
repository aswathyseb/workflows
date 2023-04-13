### Description

This folder contains ruffus variant calling workflow.


### Workflow details

The main workflow is in `vcf.py` file.

All the helper functions that the main workflow uses are in helper.py

### Workflow requirements

1. Software requirements 
   1. ruffus
   2. bwa
   3. bcftools
2. Configuration requirements
   1. All parameters are specified through a YAML configuration file, eg: config.yaml 
   2. A samplesheet in comma-separated format with atleast 2 columns with column headers as 'sample,read1,read2'


An example sample sheet would be

    sample,read1,read2
    S1,data/S1_R1.fq,data/S1_R2.fq
    S2,data/S2_R1.fq,data/S2_R2.fq`
