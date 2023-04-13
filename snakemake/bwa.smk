import pandas as pd
import os

configfile: "config.yaml"

REF_DIR = config.get('ref').get('ref_dir')
IDX_DIR = config.get('ref').get('idx_dir')
ALIGN_DIR = config.get('align').get('align_dir')
VCF_DIR = config.get('vcf').get('vcf_dir')

NCPU = config.get('NCPU')

REF = config.get('ref').get('genome')
IDX = f"{IDX_DIR}/{os.path.split(REF)[1]}"
IDX_FILE = f"{IDX}.bwt"

SAMPLES = pd.read_csv(config["sample_sheet"],sep="\t").set_index('sample',drop=False)


def get_reads_PE(wildcards):
    r1 = SAMPLES.loc[wildcards.sample, "read1"]
    r2 = SAMPLES.loc[wildcards.sample, "read2"]
    return r1, r2


def get_reads_SE(wildcards):
    return SAMPLES.loc[wildcards.sample, "read1"]


rule all_bams:
    input:
        IDX_FILE,
        expand(ALIGN_DIR + "/{sample}.bam",sample=SAMPLES.index),
        expand(ALIGN_DIR + "/{sample}.bam.bai",sample=SAMPLES.index)


rule index_genome:
    input:
        REF
    output:
        IDX_FILE
    shell:
        "bwa index -p {IDX} {input}"

if config["align_type"] == "SE":
    rule align_SE:
        input:
            IDX_FILE,
            reads=get_reads_SE
        output:
            ALIGN_DIR + "/{sample}.bam"
        threads: NCPU
        params:
            rg=r"@RG\tID:{sample}\tSM:{sample}"
        shell:
            "bwa mem -t {threads} -R '{params.rg}' {IDX} {input.reads} | samtools sort >{output}"

if config["align_type"] == "PE":
    rule align_PE:
        input:
            IDX_FILE,
            reads=get_reads_PE
        output:
            ALIGN_DIR + "/{sample}.bam"
        threads: NCPU
        params:
            rg=r"@RG\tID:{sample}\tSM:{sample}"
        shell:
            "bwa mem -t {threads} -R '{params.rg}' {IDX} {input.reads} | samtools sort >{output}"

rule index_bam:
    input:
        ALIGN_DIR + "/{sample}.bam"
    output:
        ALIGN_DIR + "/{sample}.bam.bai"
    shell:
        "samtools index {input}"
