import ruffus.cmdline as cmdline
import yaml, os
from ruffus import *


def get_args():
    """
    Get run time parameters
    """
    parser = cmdline.get_argparse(description='Call variants from data using bcftools')
    parser.add_argument('--config', default="config.yaml", help="Configuration file")
    parser.add_argument('--genome', help="Reference genome")
    parser.add_argument('--ref_dir', help="Directory that stores reference genome")
    parser.add_argument('--idx_dir', help="Directory to store genome indices")
    parser.add_argument('--aln_dir', help="Directory to store bam files")
    parser.add_argument('--vcf_dir', help="Directory to store vcf files")

    args = parser.parse_args()

    return args


def get_all_args():
    """
    # Override config file parameters if provided at runtime.
    """
    args = get_args()

    config = yaml.safe_load(open(args.config))

    args.genome = args.genome if args.genome else config.get('genome')
    args.ref_dir = args.ref_dir if args.ref_dir else config.get('ref_dir')
    args.idx_dir = args.idx_dir if args.idx_dir else config.get('idx_dir')
    args.aln_dir = args.aln_dir if args.aln_dir else config.get('aln_dir')
    args.vcf_dir = args.vcf_dir if args.vcf_dir else config.get('vcf_dir')

    # An index file
    args.idx_file = os.path.join(args.idx_dir, os.path.basename(args.genome) + ".bwt")

    # Genome index
    args.index = os.path.splitext(args.idx_file)[0]

    # sample_sheet
    args.samples = config.get('data').get('sample_sheet')

    # Alignment mode and flags
    args.lib = config.get('data').get('library')
    args.sam_flags = config.get('sam').get('sam_flags')
    args.aln_flags = config.get('aligner').get('aln_flags')

    # Variant call flags
    args.pflags = config.get('vcf').get('pileup_flags')
    args.cflags = config.get('vcf').get('call_flags')
    return args


def pipeline_dry_run():
    print("\nThe different tasks in the pipeline are \n")
    print("index_genome: Index reference genome")
    print("align_reads: Align data to reference genome")
    print("index_bam: Index bam files")
    print("call_variants: Call variants")
    print("\n")


def read_data(df, mode):
    """
    Reads data as a list of lists
    """
    store = list()
    r1, r2 = "", ""
    for i in df.index:
        r1 = df.loc[i, 'read1']
        if mode == "PE":
            r2 = df.loc[i, 'read2']
        store.append([r1, r2]) if mode == "PE" else store.append([r1])
    return store


def get_prefix_suffix(fpath, sample):
    # suffix is everything from the last underscore.
    path, name = os.path.split(fpath)
    a, ext = os.path.splitext(name)
    vals = a.split('_')
    suffix = '_' + vals[-1] + ext
    return path, suffix


def run(args):
    if args.just_print:
        pipeline_dry_run()
        pipeline_printout(verbose=3)
    else:
        cmdline.run(args, checksum_level=3)
    # pipeline_run()
