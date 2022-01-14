# viridian_simulations_workflow

A workflow to simulate SARS-CoV-2 amplicons and reads for assembly with [viridian](https://github.com/iqbal-lab-org/viridian).

# Overview

A working pipeline to generate and assemble synthetic SARS-CoV-2 amplicons. Mismatches are detected between primers and their target sequences amd if present, PCR artifacts are simulated. With probabilities of 0.5, primers in amplicons are reverted to the reference primer target sequence or amplicon sequencing coverage is reduced for the amplicon. VGsim is used to simulate a SARS-CoV-2 phylogeny from a reference sequence and neutral evolution along the tree is then simulated with phastSim. A custom script splits the simulated genomes into amplicons based on the specified artic nCoV-2019 primer scheme and the error mode is simulated. Read sequencing is simulated using art_illumina for simulated Illumina reads or Badread for Nanopore reads, at coverages determined by the previous step. Viridian_workflow is run to map the simulated reads to the reference sequence and assemble into a consensus sequence. The simulated sequences are then masked, whereby bases that are only represented by low coverage amplicons are masked with Ns.

# Installation

```
git clone https://github.com/iqbal-lab-org/viridian_simulations_workflow
cd viridian_simulations_workflow
singularity build viridian_simulations_workflow.img Singularity.def
cd ..
```

# Usage

```singularity run viridian_simulations_workflow.img snakemake --cores 1 assess_assemblies```

# Configuration

## General
* ```processes```: Number of threads to extract amplicon sequences and run ART_illumina.
* ```reference_genome```: A SARS-CoV-2 reference genome for simulations and read mapping.
* ```proportion_illumina```: Proportion of total number of simulated genomes that will become simulated illumina reads. The remaining fraction become Nanopore simulated reads.
* ```primer_scheme```: The artic nCoV-2019 primer scheme (V3, V4 or V4.1).

## VGsim:
* ```rate_file```: A file specifying per haplotype birth, death, sampling and migration rates
* ```iterations```: Maximum number of iterations to run the model.
* ```pp_population_model_file```: File specifying number, population size and contact rates of host populations to model.
* ```mg_population_model_file```: File containing matrix of contact rates between the modelled host populations.
* ```seed```: Start seed.
* ```sample_size```: Number of SARS-CoV-2 samples to draw before terminating the simulation.

## phastSim:
* ```tree_file```: Newick tree of SARS-CoV-2 simulated genomes.
* ```output_dir```: Output directory name for phastSim
* ```seed```: Start seed.

## split_amplicons:
* ```random_dropout_probability```: Probability of random dropout of amplicons to simulate PCR errors.
* ```match_coverage_mean```: Mean coverage for amplicons containing **no mismatches** in the primer-binding region.
* ```match_coverage_sd```: Standard deviation of coverage for amplicons containing **no mismatches** in the primer-binding region.
* ```mismatch_coverage_mean```: Mean coverage for amplicons containing **mismatches** in the primer-binding region.
* ```mismatch_coverage_sd```: Standard deviation of coverage for amplicons containing **mismatches** in the primer-binding region.

## simulate_reads:
* ```illumina_read_length```: Length of ART Illumina reads.

## mask_assemblies:
* ```apply_mask```: Mask low coverage bases in the simulated truth genomes.