#!/usr/bin/python3

import argparse
import numpy as np
import pandas as pd
import cooler
from math import exp
from os.path import exists
from collections import defaultdict
from scipy.signal import convolve2d
import os

parser = argparse.ArgumentParser(description = 'Make coolers from numpy matrices')
parser.add_argument('matrices', help = 'List of matrices')
parser.add_argument('bin_size', type=int, help = 'Bin sizes of matrices')
parser.add_argument('oe_vectors', help = 'Base path or file of oe vectors')
parser.add_argument('regions', help = 'Regions used for predictions')
parser.add_argument('out_path', help = 'Path to store cool files')
parser.add_argument('out', help = 'Merged cooler name')
parser.add_argument('-p', '--prefix', help = 'Prefix folder for matrices (with trailing slash)')
parser.add_argument('-c', '--convolve', action='store_true', help = 'Whether to apply convolution to smooth the matrices')
args = parser.parse_args()

if not os.path.exists(args.out_path):
    os.makedirs(args.out_path)

def fix_region_id(name):
    """
    Added this to solve finding region of the coordinate, but not entirely sure why it was necessary
    """
    chrom, start = name.split('_')
    current_diff = 1000000

    for region_id, coords in regions.items():
        if chrom == coords[0]:
            new_diff_start = abs(int(start) - coords[1])
            if new_diff_start < current_diff:
                current_diff = new_diff_start
                current_region = region_id

    return current_region

def format_counts(mat_name, mat, oe_vector):
    """
    Format the output of Cleopatra into a dataframe that cooler can use
    """
    chrom, start = mat_name.split('_')
    start = int(start)
    current_mat = np.load(mat)
    if args.convolve:
        current_mat = convolve2d(current_mat, np.ones((3, 3)) / 9, mode='same')
    out_bins = []

    # go through the matrix
    for iy, ix in np.ndindex(current_mat.shape):
        diag_idx = abs(ix - iy)
        # convert oe value predicted by Cleopatra to observed value
        oe_val = exp(current_mat[iy, ix]) - 1
        obs_val = max(oe_val*oe_vector[diag_idx], 10**(-5))
        # get start and end coordinates based on bin index
        start_coord1 = start + args.bin_size*iy
        end_coord1 = start_coord1 + args.bin_size
        start_coord2 = start + args.bin_size*ix
        end_coord2 = start_coord2 + args.bin_size
        bed2g = [chrom, start_coord1, end_coord1, chrom, start_coord2, end_coord2, obs_val]
        out_bins.append(bed2g)

    out_df = pd.DataFrame(out_bins, columns=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'count'])

    return out_df
    
def make_cooler(name, df):
    """
    Convert the dataframes to cooler
    """
    bed_df = df.merge(hg38_chrom_bins, 
             how='left', 
             left_on=['chrom1', 'start1', 'end1'],
             right_on=['chrom', 'start', 'end'])
    
    bed_df = bed_df.merge(hg38_chrom_bins, 
             how='left', 
             left_on=['chrom2', 'start2', 'end2'],
             right_on=['chrom', 'start', 'end'])
    
    bed_df_clean = bed_df.rename(columns={'bin_id_x': 'bin1_id', 'bin_id_y': 'bin2_id'})
    bed_df_clean = bed_df_clean[['bin1_id', 'bin2_id', 'count']]
    bed_df_clean.dropna(inplace=True)
    test_pixels = cooler.create.sanitize_pixels(hg38_chrom_bins, tril_action = 'drop')(bed_df_clean)

    cooler.create_cooler(name, hg38_chrom_bins, test_pixels, dtypes = {'count': 'float64'})

hg38_chrom = cooler.util.read_chromsizes('/mnt/md0/clarice/src/genomes/hg38.sorted.chrom.sizes')
hg38_chrom_bins = cooler.util.binnify(hg38_chrom, args.bin_size)
hg38_chrom_bins['bin_id'] = hg38_chrom_bins.index

regions_by_chrom = defaultdict(list)
regions = {}

with open(args.regions, 'r') as f:
    header = f.readline()
    for line in f:
        chrom, start, end, region_id = line.strip('\n').split('\t')
        if region_id.startswith('region'):
            region_id = region_id.strip('region')
        regions_by_chrom[chrom].append(region_id)
        regions[region_id] = [chrom, int(start), int(end)]

all_mats = {}
mat_region = {}

# load all matrices
with open(args.matrices, 'r') as f:
    for line in f:
        chrom, start = line.strip('\n').split('/')[-1].split('_')[1:3]
        name = '_'.join([chrom, start])
        if args.prefix != None:
            all_mats[name] = args.prefix + line.strip('\n')
        else:
            all_mats[name] = line.strip('\n')
        for region_id in regions_by_chrom[chrom]:
            coords = regions[region_id]
            if int(start) > coords[1] and int(start) < coords[2]:
                mat_region[name] = region_id
                break

cool_files = []
count = 0

# convert matrices to coolers
for name, mat in all_mats.items():
    cooler_filename = args.out_path + name + '.cool'
    cool_files.append(cooler_filename)
    count += 1
    if count%10 == 0:
        print('Number files processed: ' + str(count) + '\n')
    # only generate cooler if it doesn't already exist
    if exists(cooler_filename) == False:
        try:
            region_id = mat_region[name]
        except:
            print(name)
            region_id = fix_region_id(name)
        # different handling of OE vectors depending on whether it's Micro-C or RCMC data, Micro-C only has one file of OE vectors for the whole genome
        if os.path.exists(args.oe_vectors):
            oe_vector = np.load(args.oe_vectors)
        else:
            oe_vector = np.load(args.oe_vectors + '_' + region_id + '.npy')
        out_df = format_counts(name, mat, oe_vector)
        make_cooler(cooler_filename, out_df)

# merge all created coolers
cooler.merge_coolers(args.out, cool_files, 100000, agg={'count': 'mean'})
