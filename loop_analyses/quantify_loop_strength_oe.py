#!/usr/bin/python3

import cooler
from cooltools.api import snipping
import cooltools
import pandas as pd
import pickle
import numpy as np
import bioframe
import argparse
import warnings
warnings.simplefilter(action='ignore', category=FutureWarning)

parser = argparse.ArgumentParser(description='Quantify loops by observed/expected')
parser.add_argument('--clr', help='Cooler to quantify from')
parser.add_argument('--bin_size', type=int, help='Bin size to use')
parser.add_argument('--loops', help = 'Loops to quantify. Can be in .bedpe format or just .txt format with chrom, mid1, mid2')
parser.add_argument('--quant_size', type=int, help='Total quant size around loop')
parser.add_argument('--outfile', help='Output filename')
parser.add_argument('--regions', help='File containing regions if using RCMC')
parser.add_argument('--oe_vecs', help='Path to oe vectors if using predictions')
args = parser.parse_args()

def find_region(chrom, mid1, mid2):
    """
    Figure out which RCMC region the loop is in
    """
    if args.regions:
        for idx, region in region_df.iterrows():
            if chrom == region['chrom'] and mid1 > region['start'] and mid1 < region['end']:
                return region['name']
                break
        return 'nope'
    else:
        for idx, region in hg38_arms.iterrows():
            if chrom == region['chrom'] and mid1 > region['start'] and mid1 < region['end']:
                if mid2 > region['start'] and mid2 < region['end']:
                    return region['name']
                    break
        return 'nope'

chroms = []
with open('/mnt/md0/clarice/src/genomes/hg38.sorted.chrom.sizes', 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        if line[0] != 'chrY':
            chroms.append(line[0])

clr = cooler.Cooler(f'{args.clr}::resolutions/{args.bin_size}')

hg38_chromsizes = bioframe.fetch_chromsizes('hg38')
hg38_cens = bioframe.fetch_centromeres('hg38')
hg38_arms = bioframe.make_chromarms(hg38_chromsizes, hg38_cens)
hg38_arms = hg38_arms[hg38_arms.chrom.isin(clr.chromnames)].reset_index(drop=True)
hg38_arms = hg38_arms[hg38_arms['chrom'] != 'chrY']

if args.regions:
    region_df = pd.read_csv(args.regions, sep='\t').rename(columns={'region_id':'name'})
    region_df = bioframe.sort_bedframe(region_df, view_df=hg38_chromsizes)

# get the cell type correct since it's not consistent
if 'HCT' in args.clr:
    celltype = 'HCT-116'
elif 'GM' in args.clr:
    celltype = 'GM12878'
elif 'K562' in args.clr:
    celltype = 'K562'
elif 'H1' in args.clr:
    celltype = 'H1'

# for predicted RCMC regions
if args.oe_vecs and args.regions:

    oe_vecs = {}

    for idx, row in region_df.iterrows():
        region_id = row['name'].strip('region')
        f_oe = f'{args.oe_vecs}/{args.bin_size}/{celltype}_{region_id}.npy'
        oe_vecs[row['name']] = list(np.load(f_oe))

    expected_cis = cooltools.expected_cis(
                    clr=clr,
                    view_df=region_df,
                    smooth=True,
                    aggregate_smoothed=True,
                    smooth_sigma=0.1,
                    nproc=16, 
                    clr_weight_name=None
                )

    fake_expected_all_regions = pd.DataFrame()

    for idx, row in region_df.iterrows():
        current_df = expected_cis[expected_cis['region1'] == row['name']]
        oe = oe_vecs[row['name']]
        if args.bin_size == 500:
            pad = len(current_df) - 4000
            oe = oe + ['NaN']*pad
        elif args.bin_size == 2000:
            oe = oe[:len(current_df)]
        current_df['count.avg'] = oe        
        fake_expected_all_regions = pd.concat([fake_expected_all_regions, current_df])

    fake_expected_all_regions.reset_index(inplace=True)

    snipper = snipping.ObsExpSnipper(clr, fake_expected_all_regions, view_df=region_df, expected_value_col='count.avg', cooler_opts={'balance': False})

    region_snips = {}
    for idx, row in region_df.iterrows():
        region_snips[row['name']] = snipper.select(row['name'], row['name'])

# for predicted whole genome
elif args.oe_vecs:

    oe_vecs = {}

    for chrom in chroms:
        f_oe = open(f'{args.oe_vecs}/{celltype}/{chrom}_oe_vec_{args.bin_size}.pkl', 'rb')
        oe_vec = pickle.load(f_oe)
        oe_vecs[chrom] = [np.mean(i) for i in oe_vec]

    expected_cis = cooltools.expected_cis(
                clr=clr,
                view_df=hg38_arms,
                smooth=True,
                aggregate_smoothed=True,
                smooth_sigma=0.1,
                nproc=16,
                clr_weight_name=None
            )
        
    fake_expected_all_regions = pd.DataFrame()

    for idx, row in hg38_arms.iterrows():
        current_df = expected_cis[expected_cis['region1'] == row['name']]
        oe = oe_vecs[row['chrom']]
        pad = len(current_df) - 1000
        oe = oe + ['NaN']*pad
        current_df['count.avg'] = oe
        fake_expected_all_regions = pd.concat([fake_expected_all_regions, current_df])
        print('completed: ' + row['name'])

    fake_expected_all_regions.reset_index(inplace=True)

    snipper = snipping.ObsExpSnipper(clr, fake_expected_all_regions, view_df=hg38_arms, expected_value_col='count.avg', cooler_opts={'balance': False})
    region_snips = {}
    for idx, row in hg38_arms.iterrows():
        region_snips[row['name']] = snipper.select(row['name'], row['name'])

# for rcmc regions
elif args.regions:

    expected_cis = cooltools.expected_cis(
                    clr=clr,
                    view_df=region_df,
                    smooth=True,
                    aggregate_smoothed=True,
                    smooth_sigma=0.1,
                    nproc=16
                )
    
    snipper = snipping.ObsExpSnipper(clr, expected_cis, view_df=region_df)

    region_snips = {}
    for idx, row in region_df.iterrows():
        region_snips[row['name']] = snipper.select(row['name'], row['name'])

loop_add = int(args.quant_size/2)

outfile = open(args.outfile, 'w')
if args.regions:
    outfile.write('loop_id\tloop_strength\tregion\n')
else:
    outfile.write('loop_id\tloop_strength\n')

count = 0

with open(args.loops, 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        if len(line) == 4:
            chrom, mid1, mid2, name = line
        elif len(line) == 7:
            chrom = line[0]
            mid1 = (int(line[1]) + int(line[2]))/2
            mid2 = (int(line[4]) + int(line[5]))/2
            name = line[6]
        start1 = int(mid1) - loop_add
        end1 = int(mid1) + loop_add
        start2 = int(mid2) - loop_add
        end2 = int(mid2) + loop_add
        region = find_region(chrom, int(mid1), int(mid2)) 
        if region != 'nope':
            region_select = region_snips[region]
            try:
                snip = snipper.snip(region_select, region, region, (start1, end1, start2, end2))
            # sometimes it doesn't work and repeating it just works, not too sure why
            except:
                region_rescue = snipper.select(region, region)
                snip = snipper.snip(region_rescue, region, region, (start1, end1, start2, end2))    
            loop_strength = np.nansum(snip)
            if args.regions:
                outfile.write(f'{name}\t{loop_strength}\t{region}\n')
            else:
                outfile.write(f'{name}\t{loop_strength}\n')
            count += 1
        if count%20 == 0:
            print(str(count) + ' loops quantfied')