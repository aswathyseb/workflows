import pandas as pd

# Workflow parameters

configfile: "config.yaml"

VCF_DIR = config.get('vcf').get('vcf_dir')
SAMPLES = pd.read_csv(config.get('sample_sheet'), sep="\t").set_index('sample', drop=False)

# Load rules

include: "bwa.smk"
include: "bcftools.smk"

rule all:
    input:
        VCF_DIR + "/all.vcf.gz"

