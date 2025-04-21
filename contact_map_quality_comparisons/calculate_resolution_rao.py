#!/usr/bin/python3

import cooler
import pandas as pd
import argparse
import pickle

parser = argparse.ArgumentParser(description='Calculate resolution according to Rao et al. 2014')
parser.add_argument('--input', help='pairs file')
parser.add_argument('--regions', help='filename of regions used for RCMC')
parser.add_argument('--output_basename', help='output basename')
args = parser.parse_args()

resolutions = [50, 100, 150, 200, 250, 300, 400, 500, 800, 1000]

# Regions file is a tab-separated file with the columns chrom, start, end, region_id
# The names of the columns should be exact
region_idx = pd.read_csv(args.regions, sep = '\t')

# open file for writing
output = open(args.output_basename + '_clr_bins_proportion.txt', 'w')
output.write('bin_size\tregion_id\tbin_proportion\n')

# go through each resolution
for res in resolutions:
    # load cooler associated with the resolution
    clr = cooler.Cooler(f'/mnt/coldstorage/clarice/RCMC_analyses/merged_mcools/GM12878_merged_realigned.50.mcool::resolutions/{res}')
    # go through each region in RCMC
    for idx, row in region_idx.iterrows():
        region_id = row['region_id']
        # get the matrix of reads for each region
        mat = clr.matrix(balance=False).fetch([row['chrom'], row['start'], row['end']])
        all_bins = []
        # get all counts for each linear bin
        for row in mat:
            all_bins.append(sum(row))
        # keep only bins that have >1000 reads, consistent with the definition of resolution in Rao et al. 
        filled_counts = [i for i in all_bins if i > 1000]
        # calculate proportion of bins that have >1000 reads
        bin_proportion = len(filled_counts)/len(all_bins)
        # write output to file
        output.write(f'{res}\t{region_id}\t{bin_proportion}\n')
        print(f'{res}:{region_id} done')

output.close()
