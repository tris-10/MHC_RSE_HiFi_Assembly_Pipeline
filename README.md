# MHC RSE HiFi Assembly Pipeline

This pipeline generates _de novo_ haplotypic assemblies across the MHC using 
a combination of RSE PacBio HiFi sequences and phased microarray data.

## Overview

The pipeline is a combination of existing tools (listed below) and custom python scripts 
packaged in a Snakemake workflow.  The existing tools are all installed into conda
environments the first time the pipeline is run. The workflow is compatible with 
HPC environments and an example Slurm configuration yaml is provided.

## Data

Microarray data and the haplotypic assemblies can be found in the `assembly_data` directory.
Raw PacBio reads can be found as SRA. The EUR trio was not HiFi-level quality and was 
processed with a modified pipeline. 

## Installation

1. Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html#installing) 
1. Install [Snakemake](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html)
1. Clone the repository.
1. Create cluster configuration file and place `/home/USER/.config/snakemake/PROFILE/config.yaml`.
   See the config.yaml file in the `cluster_config` directory for an example of a Slurm setup file.
1. Open `config/pb_mhc_rse_config.yaml` and change the read cap, email and reference name settings.
1. Create, compress and index reference sequences and put them into the resource directory:
   * `resources/hg38_no_mhc.fasta.gz`: Hg38 with MHC region masked.
   * `resources/hg38_PG_1K.fasta.gz`. Hg38 with a panel of MHC sequences.  The publication uses Pangenome asssemblies, 1KG assemblies and the 7 alternative MHC haplotypes included hg38.
   * `resources/hg38_MHC_chr6.fasta`.  Hg38 chromosome 6 sequence.
   * `hg38_MHC_regions.bed`. Bed file listing the coordinates of the MHC sequences in `resources/hg38_PG_1K.fasta.gz`.
      
## How to Run

1. Copy compressed RSE HiFi fastq sequences in the `fastq` directory. The file should
   be named `SAMPLE.fastq.gz`. `SAMPLE` should be a sample-specific identifier that will used to prefix all output files.
1. Copy phased microarray data into the `microarray` directory.  Three files are expected:
    - fasta file with 60bp sequence upstream of microarray SNP. The file should be named
      `SAMPLE_front.fasta`.
    - fasta file with 60bp sequence downstream of microarray SNP. The file should be named
      `SAMPLE_back.fasta`
    - Table with microarray probe information.  The file should have five columns separated
      with tabs: SNP_ID, chromosome, position, hap1 base, hap2 base. Phasing
      can be done with Trios or predicted using a tool like ShapeIT. File should be namd
      `SAMPLE_table.txt`
1. Navigate to workflow directory and launch the pipeline: `snakemake --profile PROFILE`.
   `PROFILE` should match the directory name in step 4 of the installation.
1. The pipeline will run for several hours, so consider launching in a screen session.
1. The assemblies will be located in the `results/SAMPLE_final_contigs/` directory on successful
   completion of the pipeline. There will be one fasta file per haplotype. 
      
## Known Issues
    
* All necessary tools will be installed the first time the pipeline is launched. 
   The only exception is amos, which isn't setup properly when installed via conda.
   After the conda install, open the `minimus2` binary in `/home/USER/miniconda3/envs/CONDA_NAME/bin/minimus2` 
   and change the `SHOW_COORDS` path to ~/miniconda3/envs/CONDA_NAME/bin/show-coords. The 
   binary prefix should match the `NUCMER` prefix.

## External tools

* [minimap2](https://github.com/lh3/minimap2)
* [canu](https://github.com/marbl/canu)
* [shapeit4](https://odelaneau.github.io/shapeit4/)
* [whatshap](https://github.com/whatshap/whatshapcan)
* [longshot](https://github.com/pjedge/longshot)
* [bwa](https://github.com/lh3/bwa)
* [jellyfish](https://github.com/gmarcais/Jellyfish)
* [amos](https://amos.sourceforge.net/wiki/index.php/Minimus)
* [samtools](https://github.com/samtools/)
* [bedtools](https://github.com/arq5x/bedtools2)

