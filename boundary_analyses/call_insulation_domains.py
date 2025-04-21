#!/usr/bin/python3

import argparse
import cooler
import cooltools
import bioframe
import pandas as pd
import os

parser = argparse.ArgumentParser(description = 'Call insulation domain boundaries')
parser.add_argument('clr', help='Cooler to calculate insulation scores on')
parser.add_argument('window', help='List of window sizes to calculate on')
parser.add_argument('out_path', help='Output path')
args = parser.parse_args()

def convert_bp(bp):

    bp = int(bp)

    if bp >= 1000:
        out = str(int(bp/1000)) + 'kb'
    else:
        out = str(bp) + 'bp'

    return out

def get_compartments(df, res, outfile):
    current_boundaries = df[df['is_boundary_' + res]]
    previous_boundary = 'NA'
    previous_region = 'NA'

    with open(outfile, 'w') as f:
        for idx, row in current_boundaries.iterrows():
            current_boundary = int((row['start'] + row['end'])/2)
            current_region = row['region']
            if previous_boundary != 'NA' and current_region == previous_region:
                f.write('{chrom}\t{start}\t{end}\t{name}\t{score}\n'.format(
                    chrom=row['chrom'], start=previous_boundary+1, end=current_boundary, name='.', score=row['boundary_strength_' + res]))
                previous_boundary = current_boundary
                previous_region = current_region
            else:
                previous_boundary = current_boundary
                previous_region = current_region

if not os.path.exists(args.out_path):
    os.makedirs(args.out_path)

hg38_chromsizes = bioframe.fetch_chromsizes('hg38')
hg38_cens = bioframe.fetch_centromeres('hg38')
hg38_arms = bioframe.make_chromarms(hg38_chromsizes, hg38_cens)

region_df = pd.read_csv('/mnt/md0/clarice/src/region_idx.txt', sep='\t').rename(columns={'region_id':'name'}) 
region_df = bioframe.sort_bedframe(region_df, view_df=hg38_arms)

missing_probes = pd.read_csv('/mnt/md0/clarice/src/v1_v2_bases_not_covered.bed', sep='\t', names=['chrom', 'start', 'end'])

clr = cooler.Cooler(args.clr)

all_windows = args.window.split(',')

insulation_table = cooltools.insulation(clr, [int(window) for window in all_windows], view_df=region_df)
insulation_table = bioframe.setdiff(insulation_table, missing_probes)

for window in all_windows:
    outfile = args.out_path + 'boundaries_' + convert_bp(window) + '.bed'
    print('starting')
    get_compartments(insulation_table, window, outfile)

