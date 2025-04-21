### Scripts for loop-related analyses
Many of these will use RCMC loops as examples but the same can be applied to Cleopatra loops.

1.`calculate_loop_anchor_epi_scores.py`: calculate epigenomic enrichment on loop anchors 
The loops are a bedpe file with the size of the anchors that you want to calculate enrichment on.  
The epigenome_filenames input is a tab-separated file with two columns: 1. full path to bam file and 2. name.  

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

2. `loop_annotations.Rmd`: Annotate loops and plot loop anchor enrichments
3. `quantify_loop_strength_oe.py`: Quantify loop strength

```
usage: quantify_loop_strength_oe.py [-h] [--clr CLR] [--bin_size BIN_SIZE] [--loops LOOPS] [--quant_size QUANT_SIZE] [--outfile OUTFILE] [--regions REGIONS] [--oe_vecs OE_VECS]

Quantify loops by observed/expected

options:
  -h, --help            show this help message and exit
  --clr CLR             Cooler to quantify from
  --bin_size BIN_SIZE   Bin size to use
  --loops LOOPS         Loops to quantify. Can be in .bedpe format or just .txt format with chrom, mid1, mid2
  --quant_size QUANT_SIZE
                        Total quant size around loop
  --outfile OUTFILE     Output filename
  --regions REGIONS     File containing regions if using RCMC
  --oe_vecs OE_VECS     Path to oe vectors if using predictions
```
4. `loop_celltype_comparisons.ipynb`: 4-way comparison of loop overlaps between cell types and loop pileups
5. `loop_celltype_comparisons.Rmd`: Compare and plot loops cell-type-specific/shared loops between cell types
6. `capture_hic_comparisons.Rmd`: Compare loops between RCMC and pcHi-C from Mifsud et al. 2015
7. `capture_hic_comparisons.ipynb`: Loop pileups of RCMC/pcHi-C loops
