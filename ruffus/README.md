### Description

This folder contains ruffus variant calling workflow.


### Workflow details

The main workflow is in `vcf.py` file.

All the helper functions that the main workflow uses are in helper.py

### How to run?

    # Runs using 4 processors and parameters specified in config file.
    python vcf.py --config config.yaml  -j 4

To view all the commandline options
    
    python vcf.py -h

To perform a dry-run

    python vcf.py -n
    
To override any of the runtime parameters specified in the config file, for eg: alignment directory

    python vcf.py  --config config.yaml  -j 4 --aln_dir align

### Sample sheet specification

Sample sheet is a comma-separated file with the header sample,read1,read2

An example sample sheet would look like this:

    sample,read1,read2
    S1,data/S1_R1.fq,data/S1_R2.fq
    S2,data/S2_R1.fq,data/S2_R2.fq

### Workflow requirements

1. Software requirements 
   1. ruffus
   2. bwa
   3. bcftools
2. Configuration requirements
   1. All parameters are specified through a YAML configuration file, eg: config.yaml 
   2. A samplesheet in comma-separated format with atleast 2 columns with column headers as 'sample,read1,read2'




