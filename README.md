# viridian_simulations_workflow

A workflow to simulate SARS-CoV-2 amplicons and reads.

# Overview

A working pipeline to generate and assemble synthetic SARS-CoV-2 amplicons, with sequencing coverage conditional on the presence or absence of mismatches between primers and the target primer-binding region. VGsim is used to simulate positive selection on a SARS-CoV-2 reference sequence and neutral evolution along the tree then simulated with phastSim. A custom script splits the simulated genomes into amplicons based on the artic nCoV-2019 V3 primer scheme and sequencing coverage calculated based on whether mismatches are present between the primer and its target site. Read sequencing is simulated using art_illumina at a coverage determined by the previous step. Viridian_workflow is run to map the simulated reads to the reference sequence and assemble into a consensus sequence. A custom script then compares the viridian assemblies to the original simulated genomes output by VGsim to generate a distribution of base mismatches between the two assemblies.

# Configuration

## General
* ```processes```: Number of threads to extract amplicon sequences and run ART_illumina.
* ```reference_genome```: A SARS-CoV-2 reference genome for simulations and read mapping.
* ```proportion_illumina```: Proportion of total number of simulated genomes that will become simulated illumina reads. The remaining fraction become Nanopore simulated reads.

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
* ```primer_sequences```: TSV file of artic primer schemes with columns "name pool seq length %gc tm (use 65)"
* ```primer_positions```: BED file specifying primer start and stop positions.
* ```random_dropout_probability```: Probability of random dropout of amplicons to simulate PCR errors.
* ```match_coverage_mean```: Mean coverage for amplicons containing **no mismatches** in the primer-binding region.
* ```match_coverage_sd```: Standard deviation of coverage for amplicons containing **no mismatches** in the primer-binding region.
* ```mismatch_coverage_mean```: Mean coverage for amplicons containing **mismatches** in the primer-binding region.
* ```mismatch_coverage_sd```: Standard deviation of coverage for amplicons containing **mismatches** in the primer-binding region.

## viridian_workflow:
* ```primer_bed```: BED file specifying start and end positions of amplicons in the SARS-CoV-2 reference genome.