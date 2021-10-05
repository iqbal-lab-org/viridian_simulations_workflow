# viridian_simulations_workflow
Workflow to simulate SARS-CoV-2 amplicons and reads.

# Configuration

## General
* ```processes```: Number of threads to extract amplicon sequences and run ART_illumina.
* ```reference_genome```: A SARS-CoV-2 reference genome for simulations and read mapping.

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