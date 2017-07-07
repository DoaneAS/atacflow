#!/usr/bin/env python2

# Ashley Doane, Elemento Lab
# adapted from: Daniel Kim, CS Foo
# Script to run ataqc, all parts

import matplotlib
matplotlib.use('Agg')

import os
import sys
import pysam
import pybedtools
import metaseq
import subprocess
import multiprocessing
import timeit
import datetime
import gzip
import numpy as np
import pandas as pd
import scipy.stats
import argparse
import logging
import re
import yaml

from base64 import b64encode
from collections import namedtuple
from collections import OrderedDict
from io import BytesIO
from scipy.signal import find_peaks_cwt
from jinja2 import Template
from matplotlib import pyplot as plt
from matplotlib import mlab


# QC STUFF


# HELPER FUNCTIONS

def getFileHandle(filename, mode="r"):
    if (re.search('.gz$',filename) or re.search('.gzip',filename)):
        if (mode=="r"):
            mode="rb";
        return gzip.open(filename,mode)
    else:
        return open(filename,mode)


# QC FUNCTIONS

def determine_paired(bam_file):
    '''
    Quick function to determine if the BAM file is paired end or single end
    '''
    num_paired_reads = int(subprocess.check_output(['samtools',
                                                    'view', '-f', '0x1',
                                                    '-c', bam_file]).strip())
    if num_paired_reads > 1:
        return "Paired-ended"
    else:
        return "Single-ended"


def get_read_length(fastq_file):
    '''
    Get read length out of fastq file
    '''
    total_reads_to_consider = 1000000
    line_num = 0
    total_reads_considered = 0
    max_length = 0
    with getFileHandle(fastq_file, 'rb') as fp:
        for line in fp:
            if line_num % 4 == 1:
                if len(line.strip()) > max_length:
                    max_length = len(line.strip())
                total_reads_considered += 1
            if total_reads_considered >= total_reads_to_consider:
                break
            line_num += 1

    return int(max_length)

def get_fract_reads_in_regions(reads_bed, regions_bed):
    '''
    Function that takes in bed file of reads and bed file of regions and
    gets fraction of reads sitting in said regions
    '''
    reads_bedtool = pybedtools.BedTool(reads_bed)
    regions_bedtool = pybedtools.BedTool(regions_bed)

    reads = regions_bedtool.intersect(reads_bedtool, c=True)

    read_count = 0
    for interval in reads:
        read_count += int(interval[-1])
    fract_reads = float(read_count)/reads_bedtool.count()

    return read_count, fract_reads






def get_frip(final_bed, peaks):
    '''
    Given region sets, determine whether reads are
    falling in or outside these regions
    '''
#    logging.info('frip scores...')

    # Dnase regions
    #reads_dnase, fract_dnase = get_fract_reads_in_regions(final_bed,
    #                                                      dnase_regions)

    # Blacklist regions
    #reads_blacklist, \
    #    fract_blacklist = get_fract_reads_in_regions(final_bed,
     #                                                blacklist_regions)

    # Prom regions
    #reads_prom, fract_prom = get_fract_reads_in_regions(final_bed,
    #                                                    prom_regions)

    # Enh regions
    #reads_enh, fract_enh = get_fract_reads_in_regions(final_bed, enh_regions)

    # Peak regions
    reads_peaks, fract_peaks = get_fract_reads_in_regions(final_bed, peaks)

    return reads_peaks, \
        fract_peaks

# ===========================================================



def get_peak_counts(raw_peaks, naive_overlap_peaks=None, idr_peaks=None):
    '''
    Return a table with counts for raw peaks, IDR peaks, and naive
    overlap peaks
    '''

    # Count peaks
    raw_count = sum(1 for line in getFileHandle(raw_peaks))
    if naive_overlap_peaks != None:
        naive_count = sum(1 for line in getFileHandle(naive_overlap_peaks))
    else:
        naive_count = 0

    if idr_peaks != None:
        idr_count = sum(1 for line in getFileHandle(idr_peaks))
    else:
        idr_count = 0

    # Literally just throw these into a QC table
    results = []
    results.append(raw_count)
 #   results.append(QCGreaterThanEqualCheck('Naive overlap peaks',
 #                                          10000)(naive_count))
#    results.append(QCGreaterThanEqualCheck('IDR peaks', 10000)(idr_count))

    return raw_count



def parse_args():
    '''
    Set up the package to be run from the command line
    '''
    parser = argparse.ArgumentParser(description='ATAC-seq QC package')
    # Directories and prefixes

    parser.add_argument('--bed')#default=os.path.join(os.environ.get('${PWD}', None),os.environ.get('$Sample',None )))
    parser.add_argument('--peaks', default=None)
    #parser.add_argument('--bam', default=None)
    parser.add_argument('--out', default=None)

    args = parser.parse_args()

    #FINAL_BAM = args.finalbam
    FINAL_BED = args.bed
    PEAKS = args.peaks
    #OUTPUT_PREFIX = os.path.join(args.workdir, args.outprefix)
    #os.system('mkdir -p {0}'.format(args.outdir))
    OUTFILE = args.out
    return  FINAL_BED, PEAKS, OUTFILE

def main():

    # Parse args
    [ FINAL_BED, PEAKS, OUTFILE] = parse_args()

    reads_peaks, fract_peaks = get_frip(FINAL_BED, PEAKS)
    npeaks = get_peak_counts(PEAKS)
    #out_dir = os.path.dirname(PEAKS)
    OUTPUT_FILE = os.path.normpath(OUTFILE)
    textfile = open('{0}'.format(OUTFILE), 'w')
    textfile.write('{0}\t'.format(reads_peaks))
    textfile.write('{0}\t'.format(fract_peaks))
    textfile.write('{0}'.format(npeaks))
    textfile.close()

main()


