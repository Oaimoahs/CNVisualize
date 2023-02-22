#!/usr/bin/env python
# -*- coding: utf-8 -*-
### extract number of reads by region windows from bam files

import sys
import os
import logging
import math
import csv
import argparse
import pandas as pd


FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)

chromosomes = [str(i) for i in range(1, 23)] + ['X', 'Y']
length = [249250621, 243199373, 198022430, 191154276, 180915260, 171115067, 159138663, 146364022, 141213431, 135534747, 135006516, 133851895, 115169878, 107349540, 102531392, 90354753, 81195210, 78077248, 59128983, 63025520, 48129895, 51304566, 155270560, 59373566]


def parse_args(args):
    description = "Extract number of reads by region windows from bam files."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i","--indir", required=True, type=str, help="input bam file dir")
    parser.add_argument("-b","--barcode", required=True, type=str, help="barcode list file, each line is a unique barcode")
    parser.add_argument("-t", "--threads", type=int, required=False, default=1, help="number of threads for multi-threading")
    parser.add_argument("-o", "--outdir", type=str, required=False, default='./', help="Running directory where to write the read number csv (default: current directory)")
    args = parser.parse_args(args)

    if not os.path.isdir(args.indir):
        raise ValueError(f"Running directory does not exists: {args.indir}")
    if not os.path.isfile(args.barcode):
        raise ValueError("barcode list file does not exist!")
    # if not os.path.isdir(args.outdir):
    #     raise ValueError(f"Running directory does not exists: {args.outdir}")

    return {
        'bampath' : args.indir,
        'bc_file' : args.barcode,
        'n_threads' : args.threads,
        'rundir' : os.path.abspath(args.outdir)
    }

def get_numreads(bam, region):
    output = os.popen('samtools coverage -r %s %s' % (region, bam))
    line = output.read()
    numreads = line.strip('\n').split('\n')[1].split('\t')[3]
    return numreads

def generate_df(bam_path, out_path, barcode_list):
    for bc in barcode_list:
        bam_file = os.path.join(bam_path, bc + "_sorted.bam")
        reads = []
        logger.info('Start %s' % bc)
        for chr in chromosomes:
            logger.info('Start %s' % chr)
            reads.append('chr'+chr)
            len = length[chromosomes.index(chr)]
            for i in range(math.ceil(len/500000)):
                start = i * 500000 + 1
                end = (i + 1) * 500000
                if end > len:
                    end = len
                reg = chr + ':' + str(start) + '-' + str(end)
                numreads = get_numreads(bam_file, reg)
                reads.append(numreads)
        with open(out_path, 'a+') as f:
            csv_write = csv.writer(f)
            data_row = [bc] + reads
            csv_write.writerow(data_row)

def main(args=None):
    args = parse_args(args)
    logger.info('\n'.join(['Arguments:'] + [f'{a} : {args[a]}' for a in args]))
    
    bc_list = []
    with open(args["bc_file"]) as fin:
        for line in fin:
            barcode = line.strip()
            bc_list.append(barcode)

    generate_df(args["bampath"], args["rundir"], bc_list)

if __name__ == "__main__":
    main()