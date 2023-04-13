"""
This pipeline calls variants from the data.
"""

import csv, os, sys
from ruffus import *
import ruffus.cmdline as cmdline


def read_config(fname):
    from configparser import ConfigParser, ExtendedInterpolation
    # Create configparser object
    configs = ConfigParser(interpolation=ExtendedInterpolation())
    # Read config file
    configs.read(fname)
    samples = configs['fastq']['samples']
    suffix1 = configs['fastq']['suffix1']
    suffix2 = configs['fastq']['suffix2']
    paired = configs.getboolean('fastq', 'paired')
    aligner = configs['aligner']['tool']
    aligner_flags = configs['aligner']['aligner_flags']
    pileup_flags = configs['vcf']['pileup_flags']
    call_flags = configs['vcf']['call_flags']
    sam_flags = configs['sam']['sam_flags']
    genome = configs['runtime']['genome']
    index_dir = configs['runtime']['index_dir']
    fastq_dir = configs['runtime']['fastq_dir']
    align_dir = configs['runtime']['align_dir']
    vcf_dir = configs['runtime']['vcf_dir']
    # ncpu = int(configs['runtime']['ncpu'])

    keys = ['samples', 'suffix1', 'suffix2', 'paired', 'aligner', 'aligner_flags', 'pileup_flags', 'call_flags',
            'sam_flags', 'genome', 'index_dir', 'fastq_dir', 'align_dir', 'vcf_dir', ]
    vals = [samples, suffix1, suffix2, paired, aligner, aligner_flags, pileup_flags, call_flags, sam_flags,
            genome, index_dir, fastq_dir, align_dir, vcf_dir
            ]
    inputs = dict(zip(keys, vals))
    return inputs


def read_data(indir, fname, suffix1, suffix2):
    """
    Input is a file with sample listing , one sample per row
    Build the file name from indir, samples  and suffices.
    Returns a list of tuples with (read1, read2)
    """
    data = list()
    stream = csv.reader(open(fname), delimiter=",")

    for row in stream:
        sample = row[0]
        f1 = os.path.join(indir, sample + suffix1)
        f2 = os.path.join(indir, sample + suffix2)
        data.append([f1, f2]) if suffix2 else data.append([f1])
    return data


def pipeline_dry_run():
    print("\nThe different tasks in the pipeline are \n")
    print("index_genome: Index reference genome")
    print("align_reads: Align data to reference genome")
    print("index_bam: Index bam files")
    print("call_variants: Call variants")
    print("\n")


# Get parameters at run time

parser = cmdline.get_argparse(description='Call variants from data using bcftools')
parser.add_argument('--config', default="params.config", help="Configuration file")
parser.add_argument('--genome', help="Reference genome")
parser.add_argument('--index_dir', help="Directory to store genome indices")
parser.add_argument('--fastq_dir', help="Directory with fastq files")
parser.add_argument('--align_dir', help="Directory to store bam files")
parser.add_argument('--vcf_dir', help="Directory to store vcf files")

args = parser.parse_args()
config_file = args.config
genome = args.genome
fastq_dir = args.fastq_dir
align_dir = args.align_dir
vcf_dir = args.vcf_dir
ncpu = args.jobs
dry_run = args.just_print

# Read pipeline parameters from config file
inputs = read_config(config_file)

# Override config file parameters if provided at runtime.
genome = inputs['genome'] if not args.genome else args.genome
fastq_dir = inputs['fastq_dir'] if not args.fastq_dir else args.fastq_dir
align_dir = inputs['align_dir'] if not args.align_dir else args.align_dir
vcf_dir = inputs['vcf_dir'] if not args.vcf_dir else args.vcf_dir
idx_dir = inputs['index_dir'] if not args.index_dir else args.index_dir

# All default values from ruffus commandline
# defaults =vars(parser.parse_args([]))


# Sample list
samples = inputs['samples']

# Library mode
paired = inputs['paired']

# Suffix1 and Suffix2
s1 = inputs['suffix1']
s2 = inputs['suffix2']

# Aligner used
aligner = inputs['aligner']

# Aligner flags
aligner_flags = inputs['aligner_flags']

# Samtools flags
sam_flags = inputs['sam_flags']

# Pileup flags
pileup_flags = inputs['pileup_flags']

# Variant call flags
call_flags = inputs['call_flags']

# Get fastq files
reads = read_data(fastq_dir, samples, s1, s2)

# An index file
idx_file = os.path.join(os.path.dirname(genome), idx_dir, os.path.basename(genome) + ".bwt")

# Genome index
index = os.path.splitext(idx_file)[0]


# Task1 : Index reference genome
@transform(genome, regex(genome), idx_file)
def index_genome(inp, out):
    os.makedirs(os.path.dirname(idx_file), exist_ok=True)
    cmd = f"bwa index -p {index} {inp}"
    print(cmd)
    os.system(cmd)


# Pattern specifying read1 fastq file.
patt = f"{fastq_dir}/(.*){s1}"

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
    cmdline.run(args, checksum_level=3)
    #pipeline_run([], multiprocess=ncpu, checksum_level=3)
