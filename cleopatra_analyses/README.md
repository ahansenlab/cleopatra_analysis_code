### Scripts for Cleopatra-related analyses.

1. `make_cooler.py`: Generates coolers for holdout regions to benchmark Cleopatra. 

```
usage: make_cooler.py [-h] [-p PREFIX] [-c] matrices bin_size oe_vectors regions out_path out

Make coolers from numpy matrices

positional arguments:
  matrices              List of matrices
  bin_size              Bin sizes of matrices
  oe_vectors            Base path or file of oe vectors (RCMC has OE vectors calculated for each region, Micro-C has one OE vector for the whole genome)
  regions               Regions used for predictions
  out_path              Path to store cool files
  out                   Merged cooler name

options:
  -h, --help            show this help message and exit
  -p PREFIX, --prefix PREFIX
                        Prefix folder for matrices (with trailing slash)
  -c, --convolve        Whether to apply convolution to smooth the matrices
```

2. `make_genome_wide_cooler.py`: Make cooler from genome-wide predictions.  
Needs to be different because the predictions come out in a different format.

```
usage: make_genome_wide_cooler.py [-h] [-d DATA_PATH] [-oe OE_PATH] [-b BIN_SIZE] [-o OUT_PATH] [--mean_oe] [--use-existing]

Convert genome-wide Cleopatra predictions into cooler

options:
  -h, --help            show this help message and exit
  -d DATA_PATH, --data_path DATA_PATH
                        Path to the data
  -oe OE_PATH, --oe_path OE_PATH
                        Path to oe vectors
  -b BIN_SIZE, --bin_size BIN_SIZE
                        Cooler resolution
  -o OUT_PATH, --out_path OUT_PATH
                        Path to store cool files
  --mean_oe             Whether to take mean of inferred vectors
  --use-existing        Use existing cool files
```

3. `cleopatra_correlations.ipynb`: calculate Pearson's correlation by distance between different coolers
4. `cleopatra_correlations.Rmd`: plot Pearson's correlation by distance between different coolers
5. `less_epi_pileups.ipynb`: plot pileups for less input models
