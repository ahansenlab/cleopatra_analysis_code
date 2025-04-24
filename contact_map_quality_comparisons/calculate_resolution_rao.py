#!/usr/bin/python3

import cooler
import pandas as pd
import argparse
import pickle

parser = argparse.ArgumentParser(description='Calculate resolution according to Rao et al. 2014')
parser.add_argument('--input', help='.mcool to use for calculation')
parser.add_argument('--output_basename', help='output basename')
args = parser.parse_args()

resolutions = [50, 100, 150, 200, 250, 300, 400, 500, 800, 1000]

region_idx = pd.read_csv('/mnt/md0/clarice/src/region_idx.txt', sep = '\t')
all_mats = {}

output = open(args.output_basename + '_clr_bins_proportion.txt', 'w')
output.write('bin_size\tregion_id\tbin_proportion\n')
print('starting')

clr_name = args.input

for res in resolutions:
    clr = cooler.Cooler(f'{clr_name}::resolutions/{res}')
    all_mats[res] = {}
    for idx, row in region_idx.iterrows():
        region_id = row['region_id']
        mat = clr.matrix(balance=False).fetch([row['chrom'], row['start'], row['end']])
        # all_mats[res][row['region_id']] = mat
        all_bins = []
        for row in mat:
            all_bins.append(sum(row))
        filled_counts = [i for i in all_bins if i > 1000]
        bin_proportion = len(filled_counts)/len(all_bins)
        output.write(f'{res}\t{region_id}\t{bin_proportion}\n')
        print(f'{res}:{region_id} done')

output.close()
