#### Description
	
 Here variant calling is described using different workflows such as Ruffus, Snakemake, and GNU Make.

 All workflows use 
  
* data in `data` directory
* reference in `refs` directory
* config file and sample sheet in `config` folder.

The code for each workflow is in `src` under each workflow name.
 
For example, snakemake variant calling workflow can be invoked as
    `snakemake -s src/snakemake/variants.smk`

 Variant calling workflows included here are.
 
 1. [Ruffus Variant Caller](src/ruffus/README.md)
 2. [Snakemake Variant Caller](src/snakemake/README.md)



