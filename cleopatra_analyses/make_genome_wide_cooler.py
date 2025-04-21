import pickle
import numpy as np
import cooler
import pandas as pd
import os
import argparse

parser = argparse.ArgumentParser(description='Convert genome-wide Cleopatra predictions into cooler')
parser.add_argument('-d', '--data_path', help='Path to the data')
parser.add_argument('-oe', '--oe_path', help='Path to oe vectors')
parser.add_argument('-b', '--bin_size', type=int, help='Cooler resolution')
parser.add_argument('-o', '--out_path', help='Path to store cool files')
parser.add_argument('--mean_oe', action='store_true', help='Whether to take mean of inferred vectors')
parser.add_argument('--use-existing', action='store_true', help='Use existing cool files')
args = parser.parse_args()

if not os.path.exists(args.out_path):
    os.makedirs(args.out_path)

if not os.path.exists(args.out_path + 'coolers'):
    os.makedirs(args.out_path + 'coolers')

# find correct cell type, naming convention was not consistent
celltype = args.out_path.split('/')[-2].split('_')[0]

if 'GM' in celltype:
    oe_celltype = 'GM12878'
    celltype = 'GM12878'
elif 'HCT' in celltype:
    oe_celltype = 'HCT-116'
    celltype = 'HCT116'
else:
    oe_celltype = celltype

def format_counts(chrom, diag_idx, strata):
    """
    format Cleopatra predictions into a usable dataframe
    """
    out_bins = []
    bin_size = args.bin_size

    for idx, obs_val in enumerate(strata):
        start_coord1 = idx*bin_size
        end_coord1 = start_coord1 + bin_size
        start_coord2 = diag_idx*bin_size + idx*bin_size
        end_coord2 = start_coord2 + bin_size
        bed2g = [chrom, start_coord1, end_coord1, chrom, start_coord2, end_coord2, obs_val]
        out_bins.append(bed2g)

    out_df = pd.DataFrame(out_bins, columns=['chrom1', 'start1', 'end1', 'chrom2', 'start2', 'end2', 'count'])

    return out_df

def make_cooler(name, df):
    """
    convert dataframes into coolers
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
    test_pixels = cooler.create.sanitize_pixels(hg38_chrom_bins, tril_action = 'drop')(bed_df_clean)

    cooler.create_cooler(name, hg38_chrom_bins, test_pixels, dtypes={'count': 'float64'})

hg38_chrom = cooler.util.read_chromsizes('/mnt/md0/clarice/src/genomes/hg38.sorted.chrom.sizes')
hg38_chrom_bins = cooler.util.binnify(hg38_chrom, args.bin_size)
hg38_chrom_bins['bin_id'] = hg38_chrom_bins.index

all_chroms = list(set(hg38_chrom_bins['chrom']))
all_chroms.remove('chrY')

# Remove all normalizations
thre = np.exp(8) - 1 # The cutoff for input/output data. It's unlikely any results will be > 8, just in case.
oe_vec_thre = 1e-5
cool_files = []

for chrom in all_chroms:
    cooler_filename = args.out_path + 'coolers/' + chrom + '.cool'
    cool_files.append(cooler_filename)

    if args.use_existing and os.path.exists(cooler_filename):
        print('Using existing ' + cooler_filename)
        pass
    else:
        # Load strata for one cell type and one chromosome
        try:
            f_strata = open(f'{args.data_path}/strata/{chrom}_strata_oe.pkl', 'rb')
        except:
            f_strata = open(f'{args.data_path}/{chrom}_strata_oe.pkl', 'rb')
        strata = pickle.load(f_strata)  # This is a list of 1,000 strata

        # Load oe vectors
        f_oe = open(f'{args.oe_path}/{oe_celltype}/{chrom}_oe_vec_{args.bin_size}.pkl', 'rb')
        oe_vec = pickle.load(f_oe)

        for i in range(1000):
            # Remove log(X + 1) transform
            _s = np.where(strata[i] > thre, thre, strata[i])
            # _s = np.exp(_s) - 1
            _oe = np.where(oe_vec[i] > oe_vec_thre, oe_vec[i], oe_vec_thre)
            min_length = min(len(_s), len(_oe))
            # take the mean of all inferred OE vectors 
            if args.mean_oe:
                _s = _s[:min_length] * np.mean(_oe)
            else:
                _s = _s[:min_length] * _oe[:min_length]
            
            if i == 0:
                df = format_counts(chrom, i, _s)
            else:
                df = pd.concat([df, format_counts(chrom, i, _s)])

        make_cooler(cooler_filename, df)
        print(chrom + '_done')

# write to file
if args.bin_size == 500:
    merged_cool_name = f'{args.out_path}/{celltype}_500bp_genome_wide_predictions.cool'
elif args.bin_size == 2000:
    merged_cool_name = f'{args.out_path}/{celltype}_2kb_genome_wide_predictions.cool'

# merge all coolers
cooler.merge_coolers(merged_cool_name, cool_files, 100000, agg={'count': 'mean'})
