### Scripts for contact map quality comparisons between RCMC and Micro-C and/or Hi-C.

1. `contact_map_plotting_example.Rmd`: Example for plotting 3D genomics contact maps. 
2. `calculate_resolution_rao.py`: Calculate 'resolution' of contact maps

```
usage: calculate_resolution_rao.py [-h] [--input INPUT] [--output_basename OUTPUT_BASENAME]

Calculate resolution according to Rao et al. 2014

options:
  -h, --help            show this help message and exit
  --input INPUT         .mcool to use for calculation
  --output_basename OUTPUT_BASENAME
                        output basename
```

3. `contact_map_comparisons.ipynb`: Calculate P(s) and fraction of bins filled
4. `contact_map_comparisons_plotting.Rmd`: Plot figures. All figure plotting depends on a modified version of [plotgardener](https://phanstiellab.github.io/plotgardener/), which can be found at [this](https://github.com/claricehong/plotlandscaper) GitHub repository. 
