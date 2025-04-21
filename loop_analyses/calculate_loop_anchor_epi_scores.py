#!/usr/bin/python3

import pybedtools
import argparse

parser = argparse.ArgumentParser(description = 'Calculate the epigenome scores for each loop anchor')
parser.add_argument('bedpe', help = 'List of loops in bedpe format')
parser.add_argument('epi_files', help = 'File containing the epigenome filenames')
parser.add_argument('--outfile', '-o', help = 'Output filename')
args = parser.parse_args()

filenames = []
labels = []

# Get epigenome files
with open(args.epi_files, 'r') as f:
    for line in f:
        filename, label = line.strip('\n').split('\t')
        filenames.append(filename)
        labels.append(label)

anchors = []

# Get loop anchors
with open(args.bedpe, 'r') as f:
    for line in f:
        line = line.strip('\n').split('\t')
        anchors.append([line[0], int(line[1]), int(line[2])])
        anchors.append([line[3], int(line[4]), int(line[5])])

unique_anchors = []

# only calculate for unique anchors
for i in anchors: 
    if i not in unique_anchors:
        unique_anchors.append(i)

anchors_bt = pybedtools.BedTool(unique_anchors)
anchor_scores = anchors_bt.multicov(bams = filenames).to_dataframe(names = ['chrom', 'start', 'end'] + labels)
anchor_scores.to_csv(args.outfile, sep = '\t', index = False)
