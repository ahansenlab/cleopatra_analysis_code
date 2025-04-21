This folder contains files used for processing of RCMC data from fastq to mcool files. The pipeline generates .mcool files for each sample (in this case, cell line) with all captured regions and stats files for each replicate and capture region. Optionally, one can also generate .mcool files for each replicate, which is useful for checking concordance between replicates.  

Required packages (version number indicates ones used in this study):
    - bwa-mem2 (v2.2.1)
    - cooler (v0.10.2)
    - pairtools (v1.0.2)
    - pairix (v0.3.7)

Steps:
1. Generate list of fastq file names, then use pipeline_peripherals.ipynb to generate the sample_map.csv file for the Snakemake pipeline.
Alternatively, generate a sample_map.csv however you prefer in the following format, with one sample/rep/lane per line.
```
sample_id, sample, rep, lane, fastq_r1, fastq_r2
<sample>_<rep>_<lane>, K562, R1, L008, /path/to/fastq_R1, /path/to/fastq_R2
```

2. Gather files required for Snakemake pipeline.
    - fasta file of reference genome
    - chromsizes file of reference genome
    - .yml file of conda environment (list of packages is provided in the microc_pipeline folder, but may be easier to generate your own due to dependency issues)
    - file containing capture regions (example file provided)
    - filter_reads_merged.py (provided)

    Change the relevant paths in the Snakemake pipeline to point to your files.

3. Run the pipeline. If you also wish to generate .mcool files for each replicate, add `all_reps` to the end of the command.
```
snakemake --use-conda -s Snakefile all
```
