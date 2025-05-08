# viridian_simulations_workflow
A workflow to simulate SARS-CoV-2 amplicons and reads for assembly with viridian.

## Overview
A working pipeline to generate and assemble synthetic SARS-CoV-2 amplicons.

The main steps are outlined below:
* The pipeline begins by simulating neutral evolution from the SARS-CoV-2 reference sequence using [PhastSim](https://github.com/NicolaDM/phastSim) and a reference SARS-CoV-2 phylogeny.
* [Varifier](https://github.com/iqbal-lab-org/varifier) is then used to generate truth VCFs for each simulated assembly using the SARS-CoV-2 reference sequence.
* The primer sequences of the specified amplicon scheme are mapped to the simulated assemblies using [bwa aln](https://bio-bwa.sourceforge.net/) and the amplicon sequence extracted for each primer pair. There is then a check for mismatches between each primer and the simulated assembly. If at least one mismatch is present, PCR artifacts are simulated with equal probability: the amplicon sequencing coverage is reduced to 0 or the sequence that the primer has mapped to is reverted to the reference sequence. Random amplicon dropout is also simulated with a default probability of 0.001. The sequencing coverage of amplicons containing no artifacts is simulated using a normal distribution with mean 500 and standard deviation 20 and each amplicon is written to a FASTA file.
* Sequencing reads are simulated for each amplicon using the sequencing coverages determined in the previous step. Illumina reads are simulated using [ART illumina](https://www.niehs.nih.gov/research/resources/software/biostatistics/art/index.cfm) and Nanopore reads are simulated using [Badread](https://github.com/rrwick/Badread).
* Illumina reads are then assembled using the [Connor Lab artic nextflow pipeline](https://github.com/connor-lab/ncov2019-artic-nf) and [viridian workflow](https://github.com/iqbal-lab-org/viridian_workflow) and Nanopore reads are assembled using [Epi2me](https://github.com/epi2me-labs/wf-artic) and [viridian workflow](https://github.com/iqbal-lab-org/viridian_workflow).
* Assemblies generated with the [Connor lab pipeline](https://github.com/connor-lab/ncov2019-artic-nf) and [Epi2me](https://github.com/epi2me-labs/wf-artic) are then trimmed to remove all bases outside of the specified amplicon scheme.
* [Covid-truth-eval](https://github.com/iqbal-lab-org/covid-truth-eval) is then used to generate TSV files summarising the assembly accuracy for each tool.

## Installation
```Python
# download the repository
git clone https://github.com/iqbal-lab-org/viridian_simulations_workflow
cd viridian_simulations_workflow
# install the python dependencies in a virtual env
python3 -m venv venv
source venv/bin/activate
pip3 install -r requirements.txt
# build the containers
python3 singularity/build_images.py
```

## Inputs

The pipeline simulates assemblies from a single reference assembly that is included in this repository (`reference_genome.fasta`) along a reference SARS-CoV-2 phylogeny that is also included here (`newick_output.nwk`).

## Usage
```
snakemake --cores <threads>
```

## Configuration

### General
* ```threads```: Number of threads to extract amplicon sequences and run ART_illumina.
* ```reference_genome```: A SARS-CoV-2 reference genome for simulations and read mapping.
* ```primer_scheme```: The artic nCoV-2019 primer scheme (V3, V4 or V4.1).
* ```scheme_dir```: GitHub repository of primer schemes for the artic assembly pipeline (usually https://github.com/artic-network/primer-schemes).
* ```seed```: Start seed.
* ```container_directory```: Path to the directory containing the singularity images.
* ```nextflow_path```: Path to nextflow installation.

### phastSim:
* ```tree_file```: Newick tree of SARS-CoV-2 simulated genomes.
* ```output_dir```: Output directory name for phastSim
* ```rate_parameter```: Substitution rate for phastSim to use.

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
* ```scheme_url```: URL to the GitHub repository containing the primer scheme specification.
* ```main_nf```: Path to the Epi2me main nextflow file.
* ```illumina_workflow_container```: Path to the artic illumina singularity container.

### epi2me_badread_assemble
* ```main_nf```: Path to the Epi2me main nextflow file.
* ```nxf_sing_cache```: Directory for the nextflox cache.

### viridian_assemble:
* ```viridian_container```: Path to the viridian workflow singularity container.

## Testing
To test the functions in ```scripts/error_modes.py```, run ```python tests/run_tests.py```
