"""
This pipeline calls variants from the data.
"""

import csv, os, sys
from ruffus import *
import ruffus.cmdline as cmdline
from configparser import ConfigParser, ExtendedInterpolation


def get_args():
    """
    Get run time parameters and config file.
    """
    parser = cmdline.get_argparse(description='Call variants from data using bcftools')
    parser.add_argument('--config', default="params.config", help="Configuration file")
    parser.add_argument('--genome', help="Reference genome")
    parser.add_argument('--index_dir', help="Directory to store genome indices")
    parser.add_argument('--fastq_dir', help="Directory with fastq files")
    parser.add_argument('--align_dir', help="Directory to store bam files")
    parser.add_argument('--vcf_dir', help="Directory to store vcf files")

    args = parser.parse_args()

    # all_vals ={key: args[key] for key in vars(args)}

    config = read_config(args.config)
    config.read(args.config)



    runtime = config['runtime']

    args.genome = runtime['genome'] or args.genome
    args.idx_dir = runtime['index_dir'] or args.index_dir

    # An index file
    args.idx_file = os.path.join(os.path.dirname(args.genome), args.idx_dir, os.path.basename(args.genome) + ".bwt")
    # Genome index
    args.index = os.path.splitext(args.idx_file)[0]

    return args


def pipeline_dry_run():
    print("\nThe different tasks in the pipeline are \n")
    print("index_genome: Index reference genome")
    print("align_reads: Align data to reference genome")
    print("index_bam: Index bam files")
    print("call_variants: Call variants")
    print("\n")


def read_config(fname):
    # Create configparser object
    configs = ConfigParser(interpolation=ExtendedInterpolation())
    return configs


args = get_args()


# Task1 : Index reference genome
@transform(args.genome, regex(args.genome), args.idx_file)
def index_genome(inp, out):
    os.makedirs(os.path.dirname(out), exist_ok=True)
    index = 'refs/genome.fa'
    cmd = f"bwa index -p {index} {inp}"
    print(cmd)
    os.system(cmd)


def run():
    if args.just_print:
        pipeline_dry_run()
        pipeline_printout(verbose=3)
    else:
        cmdline.run(args, checksum_level=3)


if __name__ == "__main__":
    run()
