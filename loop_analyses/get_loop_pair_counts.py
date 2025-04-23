#!/usr/bin/python3

from collections import defaultdict
import pandas as pd
import argparse

parser = argparse.ArgumentParser(description = 'Count motif pairs in loop anchors')
parser.add_argument('celltype', help = 'Celltype')
parser.add_argument('anchor_type', help = 'Type of anchor (eg. E-P)')
args = parser.parse_args()

celltype = args.celltype
anchor_type = args.anchor_type
anchors_fimo = defaultdict(list)

with open(f'{celltype}_{anchor_type}_anchors/fimo.tsv', 'r') as f:
    header = f.readline()
    for line in f:
        line = line.strip('\n').split('\t')
        try:
            anchors_fimo[line[2]].append(line[0])
        except:
            print(line)
            continue

anchor_coords = {}

with open(f'{celltype}_{anchor_type}_anchor_coords.tsv', 'r') as f:
    header = f.readline()
    for line in f:
        line = line.strip('\n').split('\t')
        coords = '_'.join(line[1:])
        anchor_coords[coords] = line[0]

loops_with_motifs = {}
loops_coords = {}

with open(f'{celltype}_{anchor_type}_loops.tsv', 'r') as f:
    header = f.readline()
    for line in f:
        line = line.strip('\n').split('\t')
        anchor1_coords = '_'.join(line[:3])
        anchor1 = anchor_coords[anchor1_coords]
        anchor2_coords = '_'.join(line[3:6])
        anchor2 = anchor_coords[anchor2_coords]
        anchor1_motifs = anchors_fimo[anchor1]
        anchor2_motifs = anchors_fimo[anchor2]
        loops_with_motifs[line[6]] = (anchor1_motifs, anchor2_motifs)
        loops_coords[line[6]] = line[:6]

unique_shared_loops = defaultdict(dict)

with open(f'unique_shared_{anchor_type}_loops.tsv', 'r') as f:
    header = f.readline()
    for line in f:
        line = line.strip('\n').split('\t')
        unique_shared_loops[line[0]][line[1]] = line[10]

final_fimo_for_chisq = defaultdict(list)

for loop, motifs in loops_with_motifs.items():
    motif_pairs = zip(motifs[0], motifs[1])
    current_loop_motifs = []
    for motif1, motif2 in motif_pairs:
        if (motif1, motif2) not in current_loop_motifs:
            final_fimo_for_chisq[unique_shared_loops[celltype][loop]].append((motif1, motif2))
            current_loop_motifs.append((motif1, motif2))

expressed_genes = []

with open(f'{celltype}_expressed_genes.tsv', 'r') as f:
    header = f.readline()
    for line in f:
        line = line.strip('\n').split('\t')
        expressed_genes.append(line[0])

motif_pair_counts = {}

for unique_shared, loops in final_fimo_for_chisq.items():
    motif_pair_counts[unique_shared] = defaultdict(int)
    for motifs in loops:
        motif1 = motifs[0].split('.')[0]
        motif2 = motifs[1].split('.')[0]
        if motif1 in expressed_genes and motif2 in expressed_genes:
            motif_pair = [motif1, motif2]
            motif_pair.sort()
            motif_pair = tuple(motif_pair)
            if motif_pair in motif_pair_counts[unique_shared].keys():
                motif_pair_counts[unique_shared][motif_pair] += 1
            elif motif_pair in motif_pair_counts[unique_shared].keys():
                motif_pair_counts[unique_shared][motif_pair] += 1
            else:
                motif_pair_counts[unique_shared][motif_pair] += 1

with open(f'{celltype}_{anchor_type}_unique_motif_pair_counts.tsv', 'w') as f:
    f.write('motif1\tmotif2\tcount\tunique\n')
    for unique_shared, motif_counts in motif_pair_counts.items():
        for motifs, count in motif_counts.items():
            motif_1 = motifs[0]
            motif_2 = motifs[1]
            f.write(f'{motif_1}\t{motif_2}\t{count}\t{unique_shared}\n')
            