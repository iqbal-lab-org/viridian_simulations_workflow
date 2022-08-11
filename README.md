# viridian_simulations_workflow
A workflow to simulate SARS-CoV-2 amplicons and reads for assembly with viridian.

## Overview
A working pipeline to generate and assemble synthetic SARS-CoV-2 amplicons. Mismatches are detected between primers and their target sequences amd if present, PCR artifacts are simulated. With probabilities of 0.5, primers in amplicons are reverted to the reference primer target sequence or amplicon sequencing coverage is reduced for the amplicon. VGsim is used to simulate a SARS-CoV-2 phylogeny from a reference sequence and neutral evolution along the tree is then simulated with phastSim. A custom script splits the simulated genomes into amplicons based on the specified artic nCoV-2019 primer scheme and the error mode is simulated. Read sequencing is simulated using art_illumina for simulated Illumina reads or Badread for Nanopore reads, at coverages determined by the previous step. Viridian_workflow is run to map the simulated reads to the reference sequence and assemble into a consensus sequence. The simulated sequences are then masked, whereby bases that are only represented by low coverage amplicons are masked with Ns.

## Installation
```Python
# download the repository
git clone https://github.com/iqbal-lab-org/viridian_simulations_workflow
cd viridian_simulations_workflow
# install the python dependencies in a virtual env
python3 -m venv venv
source ./venv/bin/activate
pip3 install -r requirements.txt
# build the containers
python3 singularity/build_images.py
```

## Usage
```
snakemake --cores 1
```
## Configuration

### General
* ```threads```: Number of threads to extract amplicon sequences and run ART_illumina.
* ```reference_genome```: A SARS-CoV-2 reference genome for simulations and read mapping.
* ```primer_scheme```: The artic nCoV-2019 primer scheme (V3, V4 or V4.1).
* ```primer_scheme_dir```: GitHub repository of primer schemes for the artic assembly pipeline (usually https://github.com/artic-network/primer-schemes).
* ```seed```: Start seed.

### phastSim:
* ```tree_file```: Newick tree of SARS-CoV-2 simulated genomes.
* ```output_dir```: Output directory name for phastSim
* ```substitution_rate```: Substitution rate for phastSim to use.

### split_amplicons:
* ```random_dropout_probability```: Probability of random dropout of amplicons to simulate PCR errors.
* ```primer_dimer_probability```: Probability of amplicon being replaced by a primer dimer if either primer has a 3 base pair overlap with any other primer in the same PCR pool.
* ```match_coverage_mean```: Mean coverage for amplicons containing no mismatches in the primer-binding region.
* ```match_coverage_sd```: Standard deviation of coverage for amplicons containing no mismatches in the primer-binding region.
* ```mismatch_coverage_mean```: Mean coverage for amplicons containing mismatches in the primer-binding region.
* ```mismatch_coverage_sd```: Standard deviation of coverage for amplicons containing mismatches in the primer-binding region.

### simulate_reads:
* ```illumina_read_length```: Length of ART Illumina reads.

### mask_assemblies:
* ```apply_mask```: Mask low coverage bases in the simulated truth genomes.

### artic_assemble
* ```nextflow_path```: Path to nextflow installation.

### build_simulated_phylogeny, build_viridian_phylogeny, build_artic_phylogeny
* ```batch_size```: Number of samples to add to the phylogeny with UshER before running matOptimize.

## Testing
To test the functions in ```scripts/error_modes.py```, run ```python tests/run_tests.py```