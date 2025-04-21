import cooltools
import cooler
import pandas as pd
import pybedtools
import argparse

parser = argparse.ArgumentParser(description='Filter mustache loops close to blank stripes')
parser.add_argument('clr', help='cooler to find blank stripes')
parser.add_argument('mustache', help='mustache loop calls')
parser.add_argument('output', help='output filename')
parser.add_argument('-r', '--resolution', help='resolution of clr to work with')
parser.add_argument('-e', '--extend', type=int, help='number of bins to extend by')
parser.add_argument('-p', '--nproc', type=int, default=4, help='number of parallel processes')
parser.add_argument('--fdr', type=float, help='optional additional fdr threshold')
args = parser.parse_args()

pd.options.mode.chained_assignment = None

genome = '/mnt/md0/clarice/src/genomes/hg38.sorted.chrom.sizes'

clr = cooler.Cooler(args.clr + '::resolutions/' + args.resolution)
cov = cooltools.api.coverage.coverage(clr, ignore_diags=1, nproc=args.nproc, store=True)

all_chroms = ['chr' + str(i) for i in range(1, 23)]
all_chroms.append('chrX')

for idx, chrom in enumerate(all_chroms):
    bin_table = clr.bins().fetch(chrom)
    if idx == 0:
        final_bin_table = bin_table
    else:
        final_bin_table = pd.concat([final_bin_table, bin_table])

empty_mask = final_bin_table[final_bin_table['weight'].isna()]

empty_mask_bt = pybedtools.BedTool.from_dataframe(empty_mask)
empty_mask_bt = empty_mask_bt.merge()

if args.extend:
    extend_bp = int(args.resolution) * args.extend
    empty_mask_bt = empty_mask_bt.slop(l=extend_bp, r=extend_bp, g=genome)

mustache_loops = []

with open(args.mustache, 'r') as f:
    header = f.readline()
    for line in f:
        line = line.strip('\n').split('\t')
        mid1 = int((int(line[2]) + int(line[1]))/2)
        mid2 = int((int(line[5]) + int(line[4]))/2)
        chrom = line[0]
        if args.fdr:
            if float(line[6]) < args.fdr:
                mustache_loops.append(line[:6] + [f'{chrom}_{mid1}_{mid2}'])
        else:
            mustache_loops.append(line[:6] + [f'{chrom}_{mid1}_{mid2}'])
            
mustache_loops_bt = pybedtools.BedTool(mustache_loops)

mustache_loops_bt.pair_to_bed(b=empty_mask_bt, type='neither').saveas(args.output)


