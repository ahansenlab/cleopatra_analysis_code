The files in this folder are for contact map quality comparisons in Supplementary figure 2. 

Files:
1. calculate_resolution_rao.py: for calculating 'resolution' of contact maps
Usage:
```
usage: calculate_resolution_rao.py [-h] [--input INPUT] [--regions REGIONS] [--output_basename OUTPUT_BASENAME]

Calculate resolution according to Rao et al. 2014

options:
  -h, --help            show this help message and exit
  --input INPUT         pairs file
  --regions REGIONS     filename of regions used for RCMC
  --output_basename OUTPUT_BASENAME
                        output basename
```

2. contact_map_comparisons.ipynb: calculate P(s) and fraction of bins filled
3. contact_map_comparisons_plotting.Rmd: plot figures
