# HLA typing from WES data - 

An example of a simple NGS based HLA typing benchmarking study. Inspired by https://github.com/nikolasthuesen/hla-typing-benchmark

This project contains a Snakemake workflow that runs a full HLA typing pipeline, where WES samples from the 1000 Genomes project are HLA typed using Optitype, Kourami, HLA*LA and HISAT-genotype. 
The reference HLA typing is taken from DOI: 10.1371/journal.pone.0097282 
specifically, http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/working/20140725_hla_genotypes/20140702_hla_diversity.txt


This is a simplified implementation of the original benchmarking study. This project does NOT include
 - HLA typing using STC-seq
 - A study of the impact of the depth of coverage of the whole-exome sequencing sample
 - An analysis of Optitype's performance on simulated ancient DNA.
 - A detailed gold standard dataset, where newer results from doi: 10.1371/journal.pone.0206512 and doi: 10.1093/nar/gkt481 are considered

Additionally, the pipeline is set up to HLA type and evaluate the results from only two of the individuals from the 1000 Genomes dataset whereas the original study included 829 samples. Including additional samples is, however, relatively easy as the config file (snakemake/config.yaml) can relatively easily be modified. For example using create_config.py


## Software requirements
The HLA typing tools are called using Singularity which runs Docker images from https://hub.docker.com/.
No installation of the specific HLA typing tools is therefore needed.

1. [Python 3](https://www.python.org/)
2. [Singularity](https://docs.sylabs.io/guides/3.0/user-guide/installation.html)


## System requirements
- At least 40 GB of memory is needed to run the most memory-heavy step which is indexing HLA*LA's graph structure.
- The standard implementation of the Snakemake workflow currently uses 16 cores. If more/less is available, the Makefile can be easily modified.
- The WES samples are relatively large and since they are saved as both both CRAM, BAM and FASTQ files in the pipeline, a significant amount of storage is needed. However, when the results are availble for a sample, intermediate files can be deleted. Alternatively, the Snakemake script can be modified by marking large intermediate files with "temp" (see https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html#protected-and-temporary-files)


## Usage

Note, that running the pipeline will take several hours per sample, so it is recommended to do it in a tmux

1. Clone the repository

```
git clone git@github.com:nikolasthuesen/springers-hla-benchmark.git
```

2. Install the required packages

```
cd springers-hla-benchmark
make install
source virt/bin/activate
```

3. Run Snakemake workflow. 

```
make run_benchmark
```

Plots and typing results will be found in the `results` folder, which is generated when running the `Snakemake` workflow





