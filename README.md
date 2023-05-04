# HiFive Workflow

## Description
This is a `snakemake` workflow to normalize Hi-C matrices by genomic distance
using `HiFive`

## Installation 

- [ ] Make sure you have
[conda](https://docs.conda.io/en/latest/miniconda.html#latest-miniconda-installer-links)
installed.

- [ ] Then create two environments to run the pipeline:
  1. Workflow environment: 
      ```
      conda env create -n hic-workflow --file conda-workflow.lst
      ```
  2. HiFive Py2 environment:
      ```
      conda env create -n hic-hifive --file conda-hifive.lst
      ```

## Running the pipeline

- [ ] Activate the workflow environment:

  ```
  conda activate hic-workflow
  ```

- [ ] Create alias for `snakemake` so that it remains accessible even if
  environment is unloaded:

  ```
  alias snakemake=$(which snakemake)
  ```

- [ ] Edit `config.yaml` to 

  1. specify samtools binary (you can find this out by typing `which samtools`)
  2. adjust `min_interactions` that indicates the minimum number of ditags between
     two genomic regions so that they are considered to be interactiving with
     each other (you can specify multiple thresholds, which is encouraged)
  3. adjust `resolution` that determines the coarseness of the hic-map (you can
     specify here multiple resolutions which is encouraged)
  4. adjust `super_resolution_factors` that smoothes the matrix by taking into account
     the ditags of neighboring genomic regions (you can specify multiple
     factors which is encouraged)

- [ ] Create `dataset.cfg`
  1. Copy template: `cp dataset.cfg.template dataset.cfg`
  2. Fill out the settings in `[global]`, then create entries for each
     replicate/condition, each time specifying the paths for read pair 1 & 2. 

- [ ] Load HiFive environment:
  ```
  conda activate hic-hifive
  ```

- [ ] Run `snakemake`:
  ```
  snakemake -j<number of threads you can afford>
  ```

