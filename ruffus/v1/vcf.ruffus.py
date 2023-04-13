"""
This pipeline calls variants from the data.
"""

import csv, os, sys
from ruffus import *
import ruffus.cmdline as cmdline
import yaml
import pandas as pd


def get_args():
    """
    Get run time parameters and config file.
    """
    parser = cmdline.get_argparse(description='Call variants from data using bcftools')
    parser.add_argument('--config', default="config.yaml", help="Configuration file")
    parser.add_argument('--genome', help="Reference genome")
    parser.add_argument('--ref_dir', help="Directory that stores reference genome")
    parser.add_argument('--idx_dir', help="Directory to store genome indices")
    parser.add_argument('--align_dir', help="Directory to store bam files")
    parser.add_argument('--vcf_dir', help="Directory to store vcf files")

    args = parser.parse_args()

    return args


def get_all_args():
    args = get_args()
    config = yaml.safe_load(open(args.config))
    args.genome = config.get('ref').get('genome') or args.genome
    args.ref_dir = config.get('ref').get('ref_dir') or args.ref_dir
    args.idx_dir = config.get('ref').get('idx_dir') or args.idx_dir
    args.align_dir = config.get('align').get('align_dir') or args.align_dir
    args.vcf_dir = config.get('vcf').get('vcf_dir') or args.vcf_dir

    # An index file
    args.idx_file = os.path.join(args.idx_dir, os.path.basename(args.genome) + ".bwt")

    # Genome index
    args.index = os.path.splitext(args.idx_file)[0]

    # sample_sheet
    args.samples = config.get('sample_sheet')

    args.align_type = config.get('align_type')

    return args


def pipeline_dry_run():
    print("\nThe different tasks in the pipeline are \n")
    print("index_genome: Index reference genome")
    print("align_reads: Align data to reference genome")
    print("index_bam: Index bam files")
    print("call_variants: Call variants")
    print("\n")


def read_data(df, mode):
    store = list()
    r1, r2 = "", ""
    for i in df.index:
        r1 = df.loc[i, 'read1']
        if mode == "PE":
            r2 = df.loc[i, 'read2']
        store.append([r1, r2]) if mode == "PE" else store.append([r1])
    return store


args = get_all_args()

SAMPLES = pd.read_csv(args.samples, sep="\t").set_index('sample', drop=False)

reads = read_data(SAMPLES, args.align_type)


# Task1 : Index reference genome
@transform(args.genome, regex(args.genome), args.idx_file)
def index_genome(inp, out):
    os.makedirs(os.path.dirname(out), exist_ok=True)
    cmd = f"bwa index -p {args.index} {inp}"
    print(cmd)
    os.system(cmd)


# Task2 : Align data to reference genome
# @follows(index_genome)
# @transform(reads, regex(reads[0]), add_inputs(args.genome), args.align_dir + "/\\1.bam")


#
# # Pattern specifying read1 fastq file.
# patt = f"{ARGS.fastq_dir}/(.*){s1}"
#
#
# # Task2 : Align data to reference genome
# @follows(index_genome)
# @transform(reads, regex(patt), add_inputs(genome), align_dir + "/\\1.bam")
# def align_reads(inp, out):
#     os.makedirs(align_dir, exist_ok=True)
#     file1 = inp[0][0]
#     if paired:
#         # Align in paired end mode
#         file2 = inp[0][1]
#         cmd = f"bwa mem -t {ncpu} {aligner_flags} {index} {file1} {file2} | samtools view -h {sam_flags}  |samtools sort > {out}"
#     else:
#         # Align in single end mode
#         cmd = f"bwa mem -t {ncpu} {index} {file1}  | samtools sort > {out}"
#     print(cmd)
#     os.system(cmd)
#
#
# # Task3 : Index bam files.
# @follows(align_reads)
# @transform(align_reads, suffix(".bam"), ".bam.bai")
# def index_bam(inp, out):
#     cmd = f"samtools index {inp}"
#     print(cmd)
#     os.system(cmd)
#
#
# # Task4 : Call variants
# @follows(index_bam)
# @transform(align_reads, formatter(), vcf_dir + "/{basename[0]}.vcf.gz")
# def call_variants(inp, out):
#     os.makedirs(vcf_dir, exist_ok=True)
#     cmd = f"""
#              bcftools mpileup {pileup_flags} -O u -f {genome} {inp} |
#              bcftools call {call_flags} -mv -O u |
#              bcftools norm -f {genome} -d all -O u | bcftools sort -O z > {out}
#             """
#     print(cmd)
#     os.system(cmd)


if __name__ == "__main__":
    if args.just_print:
        pipeline_dry_run()
        pipeline_printout(verbose=3)
    else:
        cmdline.run(args, checksum_level=3)
