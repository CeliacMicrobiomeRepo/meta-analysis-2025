# Celiac Gut Microbiome Sequencing Meta-Analysis
Scripts supporting the 2025 [meta-analysis](https://paper.link) of celiac gut microbiome sequencing data

## About

### Publication
In this work we investigated the stool and duodenal microbiomes of individuals with celiac disease using multiple 16S rRNA sequencing datasets from across the globe. Analysis involved differential abundance analysis, comparison of diversity metrics, co-occurrence network analysis, and machine learning prediction of celiac disease. The publication can be accessed [here](https://paper.link).

### The Celiac Microbiome Repository
The Celiac Microbiome Repository (CMR) is the best effort to comprehensively combine all high throughput sequencing datasets of the gut microbiome related to celiac disease. This [publication](https://paper.link) used data from version 1.0 of The Celiac Microbiome Repository (CMR), which was up to date as of **10th September 2024**. CMR is being extended beyond this date and is accessible on the [CMR GitHub Repo](https://github.com/CeliacMicrobiomeRepo/celiac-repository/tree/main/)

### The Celiac Microbiome Repository Web App
Our [R Shiny web application](https://celiac.shinyapps.io/celiac-webapp/) draws on the CMR's data, allowing for visualisation and exploration. The code behind this site is open source: [Webapp GitHub](https://github.com/CeliacMicrobiomeRepo/celiac-webapp/tree/main)

## Contents

### 16S_preprocessing_scripts/
The `16S_preprocessing_scripts` directory contains all Python and R scripts using in the processing of sequencing data before analysis or machine learning.
 - `01_download_trim_sra.py` - Script to auto download from SRA and optionally trim adapters
 - `02_dada2_pipeline_454.R` - Script to run DADA2 (optimised for 454)
 - `02_dada2_pipeline_ion_torrent.R` - Script to run DADA2 (optimised for ion torrent)
 - `02_dada2_pipeline_paired_end.R` - Script to run DADA2 (optimised for Illumina paired end)
 - `02_dada2_pipeline_single_end.R` - Script to run DADA2 (optimised for Illumina single end)
 - `03_truncate_asvs.R` - Script to truncate ASVs to align with no overhangs
 - `04_combine_trunc_asvs_phyloseq.R` - Script to merge and filter truncated ASVs across datasets into a phyloseq object

### machine_learning_scripts/
The `machine_learning_scripts` directory contains all Python scripts used to conduct the machine learning.
 - `datasets/` - Directory with combined sets of datasets ready for ML
 - `plot_scripts/` - Directory with scripts for plotting the results of ML
 - `01_phylogseq_to_tsv.R` - Script to extracts data from a phyloseq object needed for ML
 - `02_construct_dataset_ASV.py` - Script to construct a ASV level ML dataset using output of truncate_asvs.R
 - `02_construct_dataset_taxa.py` - Script to construct a taxonomic feature ML dataset using output of DADA2 script(s)
 - `03_kfold_train.py` - Script to performs K-fold training on an ML dataset
 - `03_logo_xset_train.py` - Script to performs LOGO and XSET training on an ML dataset
 - `03_custom_sklearn.py` - Script to includes a custom MLP model that has dropout

### analysis_scripts/
The `analysis_scripts` directory contains all R scripts used to conduct the analysis.
 - `???.R` - Script to ???

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
- **Haig Bishop**:   haig.bishop@pg.canterbury.ac.nz
- **Peter Prendergast**
