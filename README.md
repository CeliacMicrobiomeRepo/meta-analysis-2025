# Celiac Gut Microbiome Sequencing Meta-Analysis
This repository contains scripts supporting the 2025 [meta-analysis](https://paper.link) of celiac gut microbiome sequencing data originating from the Celiac Microbiome Repository.

## About

### Publication
In this work we investigated the stool and duodenal microbiomes across the progession of celiac disease using multiple 16S rRNA sequencing datasets from across the globe. Analysis involved differential abundance analysis, comparison of diversity metrics, and machine learning prediction of celiac disease. The publication can be accessed [here](https://paper.link).

### The Celiac Microbiome Repository
The Celiac Microbiome Repository (CMR) is the best effort to comprehensively combine all high throughput sequencing datasets of the gut microbiome related to celiac disease. This [publication](https://paper.link) used data from version 1.0 of The Celiac Microbiome Repository (CMR), which was up to date as of **15th July 2025*. CMR is being extended beyond this date and is accessible on the [CMR GitHub Repo](https://github.com/CeliacMicrobiomeRepo/celiac-repository/tree/main/)

### The Celiac Microbiome Repository Web App
Our [R Shiny web application](https://celiac.shinyapps.io/celiac-webapp/) draws on the CMR's data, allowing for visualisation and exploration. The code behind this site is open source: [Webapp GitHub](https://github.com/CeliacMicrobiomeRepo/celiac-webapp/tree/main)

## Contents

### preprocessing/
The `preprocessing` directory contains R scripts for processing of sequencing data to prepare for analysis and machine learning. These script draw from the data in the CMR. Downloading and trimming of raw reads as well as DADA2 scripts are found in the CMR's GitHub repository.
 - `truncate_asvs.R` - Script to truncate ASVs to align with no overhangs
 - `combine_trunc_asvs_to_phyloseq.R` - Script to merge and filter truncated ASVs across datasets into phyloseq objects
 - `get_read_counts.R` - Script for summarising read counts of samples in the phyloseq objects after filtering ASVs

### analysis/
The `analysis` directory contains all R scripts used to conduct the analysis.
 - `alpha_diversity.R` - Script for alpha diversity analysis
 - `differential_abundance.R` - Script for differential abundance analysis
 - `beta_diversity.R` - Script for unconstrained beta diversity analysis
 - `constrained_beta_diversity.R` - Script for constrained beta diversity analysis (CAP)
 - `abundance_ratios.R` - Script to investigate *Firmicutes/Bacteroidota* and *Prevotella/Bacteroides* abundance ratios

### machine_learning/
The `machine_learning` directory contains all Python and R scripts used to conduct the machine learning.
 - `01_phylogseq_to_tsv.R` - Script to extract data from a phyloseq object for ML
 - `02_construct_dataset.py` - Script to construct a specific ASV-level ML dataset
 - `03_kfold_train.py` - Script to perform K-fold training on an ML dataset
 - `03_logo_xset_train.py` - Script to perform LOGO and XSET training on an ML dataset
 - `custom_sklearn.py` - Script that includes a custom MLP model with dropout

## Requirements & Licenses

The code within this repository is licensed under the **GNU Affero General Public License v3.0 (AGPL-3.0)** (see the `LICENSE` file).

The scripts rely on several external tools and libraries, each with its own license:

### Python
[Python 3.12](https://www.python.org/downloads/release/python-3120/) was used for many scripts in this repository. Packages utilised include:
 - `pandas` ([BSD 3-Clause](https://github.com/pandas-dev/pandas/blob/main/LICENSE))
 - `composition_stats` ([BSD 3-Clause](https://github.com/ntessore/composition_stats/blob/main/LICENSE.txt))
 - `numpy` ([BSD 3-Clause](https://github.com/numpy/numpy/blob/main/LICENSE.txt))
 - `sklearn` (scikit-learn) ([BSD 3-Clause](https://github.com/scikit-learn/scikit-learn/blob/main/COPYING))
 - `xgboost` ([Apache License 2.0](https://github.com/dmlc/xgboost/blob/master/LICENSE))
 - `matplotlib` ([PSF license](https://github.com/matplotlib/matplotlib/blob/main/LICENSE/LICENSE))

### R
[R version 4.4](https://www.r-project.org/) was used for many scripts in this repository. Packages utilised include:
 - `phyloseq` ([AGPL-3](https://www.bioconductor.org/packages/release/bioc/html/phyloseq.html))
 - `dada2` ([LGPL-3.0](https://github.com/benjjneb/dada2/blob/master/LICENSE))
 - `ggplot2` ([MIT](https://github.com/tidyverse/ggplot2/blob/main/LICENSE.md))
 - `Biostrings` ([Artistic License 2.0](https://bioconductor.org/packages/release/bioc/html/Biostrings.html))
 - `DECIPHER` ([GPL-3](https://bioconductor.org/packages/release/bioc/html/DECIPHER.html))

### SRAtoolkit
[SRAtoolkit](https://github.com/ncbi/sra-tools/wiki/02.-Installing-SRA-Toolkit) was used for downloading raw reads from NCBI SRA. (License: [Public Domain - U.S. Government Work](https://github.com/ncbi/sra-tools/blob/master/LICENSE))

### FastQC
[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/) was used for checking for the presence of adapter sequences. (License: [GPL v3.0](https://github.com/s-andrews/FastQC/blob/master/LICENSE.txt))

### Trimmomatic
[Trimmomatic](https://github.com/usadellab/Trimmomatic/releases) was used for trimming adapter sequences from raw reads. (License: [GPL v3.0](https://github.com/usadellab/Trimmomatic/blob/main/LICENSE))

### Cutadapt
[Cutadapt](https://cutadapt.readthedocs.io/en/stable/installation.html) was used for trimming adapter sequences from raw reads. (License: [MIT](https://github.com/marcelm/cutadapt/blob/main/LICENSE))

### DADA2
[DADA2](https://benjjneb.github.io/dada2/tutorial.html) was used for processing raw 16S reads into ASVs and abundances. (License: [LGPL-3.0](https://github.com/benjjneb/dada2/blob/master/LICENSE))

### Mothur
[Mothur](https://github.com/mothur/mothur/releases/tag/v1.48.2) was used for aligning ASVs to a full length 16S gene reference alignment. (License: [GPL v3.0](https://github.com/mothur/mothur/blob/main/LICENSE.md))


## Authors
- **Haig Bishop**:    haig.bishop@pg.canterbury.ac.nz
- **Peter Prendergast**:    peter.prendergast@pg.canterbury.ac.nz

