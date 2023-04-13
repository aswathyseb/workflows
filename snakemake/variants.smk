import pandas as pd

# Workflow parameters

configfile: "config.yaml"

VCF_DIR = config.get('vcf_dir')
SHEET = config.get('data').get("sample_sheet")
SAMPLES = pd.read_csv(SHEET,sep=",").set_index('sample',drop=False)

# Load rules

include: "bwa.smk"
include: "bcftools.smk"

rule all:
    input:
        #VCF_DIR + "/all.vcf.gz"
        expand(VCF_DIR + "/{sample}.vcf.gz",sample=SAMPLES.index),
        expand(VCF_DIR + "/{sample}.vcf.gz.csi",sample=SAMPLES.index)
