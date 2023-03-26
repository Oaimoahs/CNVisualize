#!/usr/bin/env python
# -*- coding: utf-8 -*-
### extract number of reads by region windows from bam files

import os
import logging
import csv
import argparse
from sklearn import preprocessing


FORMAT = '%(levelname)s %(asctime)-15s %(name)-20s %(message)s'
logging.basicConfig(level=logging.INFO, format=FORMAT)
logger = logging.getLogger(__name__)

def parse_args(args):
    description = "Extract number of reads by region windows from bam files."
    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("-i","--indir", required=True, type=str, help="input bam file dir")
    parser.add_argument("-b","--barcode", required=True, type=str, help="barcode list file, each line is a unique barcode")
    parser.add_argument("-c","--chromsizes", required=True, type=str, help="chromosome size file for the reference genome, each line is a chromosome")
    parser.add_argument("-w", "--windowsize", type=int, required=False, default=500000, help="fixed window size for binning the genome")
    parser.add_argument("-o", "--outdir", type=str, required=False, default='./', help="Running directory where to write the read number csv (default: current directory)")
    args = parser.parse_args(args)

    if not os.path.isdir(args.indir):
        raise ValueError(f"Running directory does not exists: {args.indir}")
    if not os.path.isfile(args.barcode):
        raise ValueError("barcode list file does not exist!")
    if not os.path.isfile(args.chromsizes):
        raise ValueError("chromosome size file does not exist!")

    return {
        'bampath' : args.indir,
        'bc_file' : args.barcode,
        'chrom_size' : args.chromsizes,
        'window_size' : args.windowsize,
        'rundir' : os.path.abspath(args.outdir)
    }

def read_chrom_size(path):
    chroms = open(path, "r")
    chrom_size_dict = {}
    for line in chroms.readlines():
        line = line.strip().split("\t")
        chrom_size_dict[line[0][3:]] = int(line[1])
    return chrom_size_dict

def get_numreads(bam, region):
    output = os.popen('samtools coverage -r %s %s' % (region, bam))
    line = output.read()
    numreads = line.strip('\n').split('\n')[1].split('\t')[3]
    return numreads

def window_indices(size, chrom_sizes):
    chroms = list(chrom_sizes.keys())
    windows_region = []
    for chrom in chroms:
        pos = 0
        while pos + 1 + size < chrom_sizes[chrom]:
            windows_region.append(chrom + ":" + str(pos+1) + "-" + str(pos+size))
            pos += size
        windows_region.append(chrom + ":" + str(pos+1) + "-" + str(chrom_sizes[chrom]))
    return windows_region

def generate_df(bc_path, bam_path, regions, out_path):
    for cell in open(bc_path, "r").readlines():
        logger.info('Start %s' % cell.strip())
        count = []
        path = os.path.join(bam_path, cell.strip() + "_sorted.bam")
        for region in regions:
            numreads = get_numreads(path, region)
            count.append(numreads)
        with open(out_path, 'a+') as f:
            csv_write = csv.writer(f)
            data_row = [cell.strip()] + count
            csv_write.writerow(data_row)

def norm_num_reads(df_raw, regions_dict, save = False):
    """
    This function normalizes the raw read counts by chromosomes

    Parameters:
        df_raw (dataframe): a dataframe of the raw read counts
        regions_dict (dictionary): a dictionary of the window numbers for each chromosome
        save (boolean): whetehr to save the normalized dataframe to *.csv or not
    """
    min_max_scaler = preprocessing.MinMaxScaler() # Min-Max normalization 
    df_norm = df_raw
    for index_cell, cell in df_raw.iterrows():
        position = 1 # skipping the first (the first is the barcode)
        for chrom in regions_dict.keys():
            # dividing read counts by chromosomes
            reads = list(cell[position : position + regions_dict[chrom]])
            # normalizing read counts
            norm_reads = min_max_scaler.fit_transform(np.array(reads).reshape(-1, 1))
            norm_reads = list(norm_reads.flatten())
            # writing into dataframe
            df_norm.loc[index_cell, position : position + regions_dict[chrom] - 1] = norm_reads
            position += regions_dict[chrom]
    if save == True:
        df_norm.to_csv('df_norm.csv', index = False)
    return df_norm

def main(args=None):
    args = parse_args(args)
    logger.info('\n'.join(['Arguments:'] + [f'{a} : {args[a]}' for a in args]))
    
    chrom_size_dict = read_chrom_size(args["chrom_size"])
    windows_region = window_indices(args["window_size"], chrom_size_dict)
    generate_df(args["bc_file"], args["bampath"], windows_region, args["rundir"])
    logger.info('Dataframe generated in %s' % args["rundir"])

if __name__ == "__main__":
    main()