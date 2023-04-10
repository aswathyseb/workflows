"""
This pipeline calls variants from the data.
"""

import csv, os, sys
from ruffus import *
import ruffus.cmdline as cmdline
from configparser import ConfigParser, ExtendedInterpolation
import helper

ARGS = helper.get_args()


# Task1 : Index reference genome
@transform(ARGS.genome, regex(ARGS.genome), ARGS.idx_file)
def index_genome(inp, out):
    print(ARGS.idx_file)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    index = 'refs/genome.fa'
    cmd = f"bwa index -p {index} {inp}"
    print(cmd)
    os.system(cmd)


if __name__ == "__main__":
    helper.run(ARGS)
