### Scripts for identifying and analyzing fine-scale insulation boundaries in 3D genomics data.

Files:
1. `call_insulation_domains.py`: get insulation boundary bed files and insulation scores

```
usage: call_insulation_domains.py [-h] clr window out_path

Call insulation domain boundaries

positional arguments:
  clr         Cooler to calculate insulation scores on
  window      List of window sizes to calculate on
  out_path    Output path

options:
  -h, --help  show this help message and exit
```

An example regions file can be found in the `src/` folder.

To convert the scores to bigwig format for plotting:
```
bedtools sort -header -i <scores_output.bedGraph> > <scores_output_sorted.bedGraph>

bedGraphToBigWig <scores_output_sorted.bedGraph> <chrom.sizes file> <scores.bw>
```

2. `boundary_analyses.Rmd`: Perform boundary-related analyses in manuscript
3. `plot_boundaries_epigenomics.Rmd`: ChIP-seq metaplots of boundaries
