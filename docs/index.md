<div align="center">

  <h1 style="font-size: 250%">istari 🔬</h1>

  <b><i>long pipeline name</i></b><br>
  <a href="https://github.com/OpenOmics/istari/actions/workflows/main.yaml">
    <img alt="tests" src="https://github.com/OpenOmics/istari/workflows/tests/badge.svg">
  </a>
  <a href="https://github.com/OpenOmics/istari/actions/workflows/docs.yml">
    <img alt="docs" src="https://github.com/OpenOmics/istari/workflows/docs/badge.svg">
  </a>
  <a href="https://github.com/OpenOmics/istari/issues">
    <img alt="GitHub issues" src="https://img.shields.io/github/issues/OpenOmics/istari?color=brightgreen">
  </a>
  <a href="https://github.com/OpenOmics/istari/blob/main/LICENSE">
    <img alt="GitHub license" src="https://img.shields.io/github/license/OpenOmics/istari">
  </a>

  <p>
    This is the home of the pipeline, istari. Its long-term goals: to accurately ...insert goal, to infer ..insert goal, and to boldly ..insert goal like no pipeline before!
  </p>

</div>  


## Overview
Welcome to istari's documentation! This guide is the main source of documentation for users that are getting started with the [whole genome regression modeling pipeline](https://github.com/OpenOmics/istari/).

The **`./istari`** pipeline is composed several inter-related sub commands to setup and run the pipeline across different systems. Each of the available sub commands perform different functions:

 * [<code>istari <b>run</b></code>](usage/run.md): Run the istari pipeline with your input files.
 * [<code>istari <b>unlock</b></code>](usage/unlock.md): Unlocks a previous runs output directory.
 * [<code>istari <b>cache</b></code>](usage/cache.md): Cache remote resources locally, coming soon!

**istari** is comprehensive pipeline that performs error-prone data preparation steps for genome-wide association studies (GWAS) optimized for WES and WGS. It relies on technologies like [Singularity<sup>1</sup>](https://singularity.lbl.gov/) to maintain the highest-level of reproducibility. The pipeline consists of a series of data processing and quality-control steps orchestrated by [Snakemake<sup>2</sup>](https://snakemake.readthedocs.io/en/stable/), a flexible and scalable workflow management system, to submit jobs to a cluster.

The pipeline is compatible with data generated from Illumina short-read sequencing technologies. As input, it accepts a VCF file (must have corresponding index file), phenotype file, and covariate file and can be run locally on a compute instance or on-premise using a cluster. A user can define the method or mode of execution. The pipeline can submit jobs to a cluster using a job scheduler like SLURM (more coming soon!). A hybrid approach ensures the pipeline is accessible to all users.

Before getting started, we highly recommend reading through the [usage](usage/run.md) section of each available sub command.

For more information about issues or trouble-shooting a problem, please checkout our [FAQ](faq/questions.md) prior to [opening an issue on Github](https://github.com/OpenOmics/istari/issues).

## Contribute

This site is a living document, created for and by members like you. istari is maintained by the members of NCBR and is improved by continous feedback! We encourage you to contribute new content and make improvements to existing content via pull request to our [GitHub repository :octicons-heart-fill-24:{ .heart }](https://github.com/OpenOmics/istari).


## References
<sup>**1.**  Kurtzer GM, Sochat V, Bauer MW (2017). Singularity: Scientific containers for mobility of compute. PLoS ONE 12(5): e0177459.</sup>  
<sup>**2.**  Koster, J. and S. Rahmann (2018). "Snakemake-a scalable bioinformatics workflow engine." Bioinformatics 34(20): 3600.</sup>  
