Scripts for loop-related analyses. Many of these will use RCMC loops as examples but the same can be applied to Cleopatra loops.

1. calculate_loop_anchor_epi_scores.py: calculate epigenomic enrichment on loop anchors 
The loops are a bedpe file with the size of the anchors that you want to calculate enrichment on.  
The epigenome_filenames input is a tab-separated file with two columns: 1. full path to bam file and 2. name.  

Usage:
```
usage: calculate_loop_anchor_epi_scores.py [-h] [--outfile OUTFILE] bedpe epi_files

Calculate the epigenome scores for each loop anchor

positional arguments:
  bedpe                 List of loops in bedpe format
  epi_files             File containing the epigenome filenames

options:
  -h, --help            show this help message and exit
  --outfile OUTFILE, -o OUTFILE
                        Output filename
```

2. loop_annotations.Rmd: for annotating loops and plotting loop anchor enrichments
3. loop_celltype_comparisons.ipynb: 4-way comparison of loop overlaps between cell types and loop pileups
4. loop_celltype_comparisons.Rmd: Compare and plot loops cell-type-specific/shared loops between cell types
5. capture_hic_comparisons.Rmd: Compare loops between RCMC and pcHi-C from Mifsud et al. 2015
6. capture_hic_comparisons.ipynb: Loop pileups of RCMC/pcHi-C loops
