import pandas as pd


configfile: "config.yaml"

ALIGN_DIR = config.get('align')['align_dir']
VCF_DIR = config.get('vcf').get('vcf_dir')
REF = config.get('ref').get('genome')

SAMPLES =pd.read_csv(config["sample_sheet"], sep="\t").set_index('sample', drop=False)

if config.get("vcf").get("call_type") == "multi-sample":
    rule call_variants:
        input:
            VCF_DIR + "/all.vcf.gz"
if config.get("vcf").get("call_type") == "sample-vcf":
    rule call_variants:
        input:
            expand(VCF_DIR + "/{sample}.vcf.gz", sample=SAMPLES.index)


rule call_multi_sample_variants:
    input:
        fa = REF,
        bam = expand(ALIGN_DIR + "/{sample}.bam", sample=SAMPLES.index),
        bai = expand(ALIGN_DIR + "/{sample}.bam.bai", sample=SAMPLES.index)
    output:
        VCF_DIR + "/all.vcf.gz"
    params:
        pflags = "-d 100",
        cflags = "--ploidy 2"
    shell:
        "bcftools mpileup {params.pflags} -O u -f {input.fa} {input.bam} | "
        "bcftools call {params.cflags} -mv -O u | "
        "bcftools norm -f {input.fa} -d all -O u | bcftools sort -O z > {output}"


rule call_sample_variants:
    input:
        fa = REF,
        bam = ALIGN_DIR + "/{sample}.bam",
        bai = ALIGN_DIR + "/{sample}.bam.bai"
    output:
        VCF_DIR + "/{sample}.vcf.gz"
    params:
        pflags = config.get('vcf')['pileup_flags'],
        cflags = config.get('vcf')['call_flags']
    shell:
        "bcftools mpileup {params.pflags} -O u -f {input.fa} {input.bam} | "
        "bcftools call {params.cflags} -mv -O u | "
        "bcftools norm -f {input.fa} -d all -O u | bcftools sort -O z > {output}"


