#!/usr/bin/python3

import argparse
import gzip as gz

parser = argparse.ArgumentParser(description = 'Filter pairs files for captured regions')
parser.add_argument('input', help = 'input pairs file')
parser.add_argument('locs', help = 'file containing capture locations')
parser.add_argument('-o', help = 'outfile directory')
parser.add_argument('-e', action = 'store_true', help = 'filter by pairs with at least one end in region')
args = parser.parse_args()

# basename for the output filename
outfile_basename = args.input.strip('.pairs.gz')

all_locs = []
open_files = {}

# read in regions used in RCMC capture
# opens one file for each region
with open(args.locs, 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        all_locs.append([line[0], line[1], int(line[2]), int(line[3])])
        current_outfile = outfile_basename.split('/')[-1] + '_loc' + line[0] + '.pairs'
        open_files[line[0]] = open(args.o + current_outfile, 'w')

def check_read(coord1, coord2, start, end):
    """
    Determines if a pair of coordinates falls within a specified capture region.

    This function checks whether both ends of a read (represented by `coord1` and `coord2`)
    are within the region defined by `start` and `end`. The function also distinguishes
    between the orientation of the pair.

    Parameters
    ----------
    coord1 : int
        The first coordinate from the pairs file.
    coord2 : int
        The second coordinate from the pairs file.
    start : int
        The start coordinate of the region to be checked (inclusive).
    end : int
        The end coordinate of the region to be checked (inclusive).

    Returns
    -------
    str or None
        Returns 'good' if both coordinates are within the region and coord1 < coord2.
        Returns 'flip' if both coordinates are within the region and coord1 > coord2.
        Returns None if the pair does not fall within the region.
    """
    if coord1 < coord2:
        if coord1 >= start and coord2 <= end:
            return 'good'
    elif coord1 > coord2:
        if coord2 >= start and coord1 <= end:
            return 'flip'
        
def check_read_one_end(coord1, coord2, start, end):
    """
    Determines if a pair of coordinates falls within a specified capture region when filtering for only one end.

    This function checks at least one end of a read (represented by `coord1` and `coord2`)
    are within the region defined by `start` and `end`. The function also distinguishes
    between the orientation of the pair.

    Parameters
    ----------
    coord1 : int
        The first coordinate from the pairs file.
    coord2 : int
        The second coordinate from the pairs file.
    start : int
        The start coordinate of the region to be checked (inclusive).
    end : int
        The end coordinate of the region to be checked (inclusive).

    Returns
    -------
    str or None
        Returns 'good' if at least one coordinate is within the region and coord1 < coord2.
        Returns 'flip' if at least one coordinate is within the region and coord1 > coord2.
        Returns None if the pair does not fall within the region.
    """
    if coord1 < coord2:
        if coord1 >= start and coord1 <= end: 
            return 'good'
        elif coord2 >= start and coord2 <= end:
            return 'good'
    elif coord1 > coord2:
        if coord1 >= start and coord1 <= end: 
            return 'flip'
        elif coord2 >= start and coord2 <= end:
            return 'flip'

# Open input pairs file
with gz.open(args.input, 'rt') as f:
    for line in f:
        # Skip header lines of pairs file
        if line.startswith('#'):
            if line.startswith('#chromsize:') or line.startswith('#samheader: @SQ'):
                pass
            else:
                # Only write relevant header lines to file
                for outfile in open_files.values():
                    outfile.write(line)
        else:
            broken_line = line.strip('\n').split('\t')
            # go through all regions
            for loc in all_locs:
                if broken_line[1] == loc[1] and broken_line[3] == loc[1]:
                    if not args.e:
                        outcome = check_read(int(broken_line[2]), int(broken_line[4]), loc[2], loc[3])
                        # write line to file if coordinates are in the region so that the first coordinate is always smaller than the second in the pair
                        if outcome == 'good':
                            open_files[loc[0]].write(line)
                        elif outcome == 'flip':
                            flipped_list = broken_line[:2] + [broken_line[4], broken_line[3], broken_line[2]] + broken_line[5:]
                            open_files[loc[0]].write('\t'.join(flipped_list) + '\n')
                    elif args.e:
                        outcome = check_read_one_end(int(broken_line[2]), int(broken_line[4]), loc[2], loc[3])
                        if outcome == 'good':
                            open_files[loc[0]].write(line)
                        elif outcome == 'flip':
                            flipped_list = broken_line[:2] + [broken_line[4], broken_line[3], broken_line[2]] + broken_line[5:]
                            open_files[loc[0]].write('\t'.join(flipped_list) + '\n')
