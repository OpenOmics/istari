<div align="center">

  <h1>istari 🔬</h1>

  **_An awesome whole genome regression modeling pipeline_**

  [![tests](https://github.com/OpenOmics/istari/workflows/tests/badge.svg)](https://github.com/OpenOmics/istari/actions/workflows/main.yaml) [![docs](https://github.com/OpenOmics/istari/workflows/docs/badge.svg)](https://github.com/OpenOmics/istari/actions/workflows/docs.yml) [![GitHub issues](https://img.shields.io/github/issues/OpenOmics/istari?color=brightgreen)](https://github.com/OpenOmics/istari/issues)  [![GitHub license](https://img.shields.io/github/license/OpenOmics/istari)](https://github.com/OpenOmics/istari/blob/main/LICENSE)

  <i>
    This is the home of the pipeline, istari. Its long-term goals: to perform error-prone data preparation for whole genome regression modeling using regenie like no pipeline before!
  </i>
</div>

## Overview
Welcome to istari! Before getting started, we highly recommend reading through [istari's documentation](https://openomics.github.io/istari/).

The **`./istari`** pipeline is composed several inter-related sub commands to setup and run the pipeline across different systems. Each of the available sub commands perform different functions:

 * [<code>istari <b>run</b></code>](https://openomics.github.io/istari/usage/run/): Run the istari pipeline with your input files.
 * [<code>istari <b>unlock</b></code>](https://openomics.github.io/istari/usage/unlock/): Unlocks a previous runs output directory.
 * [<code>istari <b>cache</b></code>](https://openomics.github.io/istari/usage/cache/): Cache remote resources locally, coming soon!

**istari** is a comprehensive pipeline that performs error-prone data preparation steps for genome-wide association studies (GWAS) optimized for WES and WGS. It relies on technologies like [Singularity<sup>1</sup>](https://singularity.lbl.gov/) to maintain the highest-level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>2</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster.

The pipeline is compatible with data generated from Illumina short-read sequencing technologies. This pipeline prepares  data for GWAS by converting file formats, performing QC filtering, preparing phenotypes, covariates, and files for [regenie](https://rgcgithub.github.io/regenie/). As input, it accepts a gVCF file (note: file must have corresponding index file), phenotype file, and covariates file. The pipeline will validate phenotype and covariates files (e.g. check file format, and add ssex and ancestry information to the covariates file if not included). Then, the pipeline performs pre-association test QC filtering for regenie step 1 based on minor allele frequency (MAF), minor allele count (MAC), linkage disequilibrium (LD), and genotype and sample missingness. The regenie developers do not recommend using >1M SNPs for step 1.  By default, certain QC thresholds are set, but QC thresholds and pruning settings can be adjusted.  The pipeline creates the annotation file, mask file, and list file, which are used as input for gene-based tests in regenie step 2.  The pipeline can submit jobs to a cluster using a job scheduler like SLURM (more coming soon!). A hybrid approach ensures the pipeline is accessible to all users.

Before getting started, we highly recommend reading through the [usage](https://openomics.github.io/istari/usage/run/) section of each available sub command.

For more information about issues or trouble-shooting a problem, please checkout our [FAQ](https://openomics.github.io/istari/faq/questions/) prior to [opening an issue on Github](https://github.com/OpenOmics/istari/issues).

## Dependencies
**Requires:** `singularity=3.10.5`  `snakemake=7.19.1`

At the current moment, the pipeline uses a mixture of environment modules and docker images; however, this will be changing soon! In the very near future, the pipeline will only use docker images. With that being said, [snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) and [singularity](https://singularity.lbl.gov/all-releases) must be installed on the target system. Snakemake orchestrates the execution of each step in the pipeline. To guarantee the highest level of reproducibility, each step of the pipeline will rely on versioned images from [DockerHub](https://hub.docker.com/orgs/nciccbr/repositories). Snakemake uses singularity to pull these images onto the local filesystem prior to job execution, and as so, snakemake and singularity will be the only two dependencies in the future.

## Installation
Please clone this repository to your local filesystem using the following command:
```bash
# Clone Repository from Github
git clone https://github.com/OpenOmics/istari.git
# Change your working directory
cd istari/
# Add dependencies to $PATH
# Biowulf users should run
module load snakemake singularity
# Get usage information
./istari -h
```

## Contribute
This site is a living document, created for and by members like you. istari is maintained by the members of OpenOmics and is improved by continuous feedback! We encourage you to contribute new content and make improvements to existing content via pull request to our [GitHub repository](https://github.com/OpenOmics/istari).

## References
<sup>**1.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**2.**  Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600.</sup>  
