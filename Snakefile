from collections import Counter
import glob
from joblib import Parallel, delayed
import math
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import os
import pandas as pd
import pickle
import random
import subprocess
import sys
from tqdm import tqdm

from primer_snps import get_primer_positions, split_amplicons, determine_coverage, find_primer_positions
from compare_amplicons import trim_genome, count_errors

configfile: 'config.yml'

rule VGsim_tree:
    input:
        rate_file=config['VGsim']['rate_file'],
        pp_population_model_file=config['VGsim']['pp_population_model_file'],
        mg_population_model_file=config['VGsim']['mg_population_model_file']
    output:
        'newick_output.nwk'
    params:
        seed=config['VGsim']['seed'],
        iterations=config['VGsim']['iterations'],
        sample_size=config['VGsim']['sample_size']
    shell:
        'python VGsim/vgsim.py {input.rate_file} -it {params.iterations} -pm {input.pp_population_model_file} \
                {input.mg_population_model_file} -seed {params.seed} -nwk --sampleSize {params.sample_size}'

rule phastSim_evolution:
    input:
        tree_file=rules.VGsim_tree.output,
        reference_genome=config["reference_genome"]
    output:
        directory(config['phastSim']['output_dir'])
    params:
        seed=config['phastSim']['seed']
    shell:
        'mkdir {output} && phastSim --outpath {output}/ --seed {params.seed} --createFasta \
                --createInfo --createNewick --createPhylip --treeFile {input.tree_file} --scale 0.333333333 \
                --invariable 0.1 --alpha 0.01 --omegaAlpha 0.01 --hyperMutProbs 0.01 0.01 --hyperMutRates 2.0 20.0 --codon \
                --scale 0.333333 --reference {input.reference_genome}'

rule split_sequences:
    input:
        rules.phastSim_evolution.output
    output:
        directory('simulated_genomes')
    run:
        with open(os.path.join(input[0], 'sars-cov-2_simulation_output.fasta')) as seqFile:
            genomes = seqFile.read().split('>')[1:]
        if not os.path.exists(output[0]):
            os.mkdir(output[0])
        for sequence in genomes:
            with open(os.path.join(output[0], sequence.splitlines()[0] + '.fasta'), 'w') as outSeq:
                outSeq.write('>' + sequence)

rule split_amplicons:
    input:
        rules.split_sequences.output
    output:
        directory('amplicon_sequences')
    params:
        primer_sequences=config['split_amplicons']['primer_sequences'],
        primer_positions=config['split_amplicons']['primer_positions'],
        random_dropout_probability=config['split_amplicons']['random_dropout_probability'],
        match_coverage_mean=config['split_amplicons']['match_coverage_mean'],
        match_coverage_sd=config['split_amplicons']['match_coverage_sd'],
        mismatch_coverage_mean=config['split_amplicons']['mismatch_coverage_mean'],
        mismatch_coverage_sd=config['split_amplicons']['mismatch_coverage_sd'],
        mismatch_probability=config['split_amplicons']['mismatch_probability'],
        processes=config['processes']
    run:
        primer_df = pd.read_csv(params.primer_sequences, sep='\t')
        # import bed file as list
        with open(params.primer_positions, "r") as pos:
            primer_metadata = pos.read().splitlines()
        # combine tsv and bed file data into dfs of left and right primers and returns list of primer sequence pairs
        left_primers, right_primers, truth_seqs = get_primer_positions(primer_metadata, primer_df)
        # create output directory
        output_dir = output[0]
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)
        # split input genomic sequences into amplicons by left primer start and right primer end
        sequence_files = glob.glob(os.path.join(input[0], '*.fasta'))
        sys.stderr.write('\nSplitting simulated sequences into amplicons\n')
        observed_primer_targets = {}
        # parallelise amplicon splitting
        job_list = [
            sequence_files[i:i + params.processes] for i in range(0, len(sequence_files), params.processes)
        ]
        for subset in tqdm(job_list):
            sample_primer_targets = Parallel(n_jobs=params.processes, prefer="threads")(delayed(split_amplicons)(genomic_sequence,
                                                                                                                left_primers,
                                                                                                                right_primers,
                                                                                                                params.mismatch_probability,
                                                                                                                output_dir) for genomic_sequence in subset)
            for elem in sample_primer_targets:
                observed_primer_targets.update(elem)
        # detect if there are any snps in the primer binding regions in our simulated samples
        sys.stderr.write("\nCalculating amplicon coverages\n")
        sample_coverages = {}
        for sample in tqdm(observed_primer_targets):
            primer_snps = {}
            for pair in truth_seqs:
                if not observed_primer_targets[sample][pair] == truth_seqs[pair]:
                    primer_snps[pair] = 1
                else:
                    primer_snps[pair] = 0
            # generate dictionary of amplicon coverages, considering mismatches and random amplicon dropout
            amplicon_coverages = determine_coverage(primer_snps,
                                                    params.random_dropout_probability,
                                                    params.match_coverage_mean,
                                                    params.match_coverage_sd,
                                                    params.mismatch_coverage_mean,
                                                    params.mismatch_coverage_sd)
            sample_coverages[sample] = amplicon_coverages
        # Pickle the output
        sys.stderr.write("\nPickling the output\n")
        with open(os.path.join(output_dir, 'amplicon_coverages.pickle'), 'wb') as ampOut:
            # Pickle using the highest protocol available.
            pickle.dump(sample_coverages, ampOut, pickle.HIGHEST_PROTOCOL)

rule simulate_reads:
    input:
        rules.split_amplicons.output
    output:
        directory('read_output')
    params:
        processes=config['processes'],
        prop_illumina=config['proportion_illumina']
    run:
        def simulate_ART_reads(genome,
                               output,
                               sample_coverages):
            """Function to run ART on amplicon sequences per simulated genomic sequence"""
            sample_name = os.path.basename(genome[:-1])
            output_dir = os.path.join(output, sample_name)
            if not os.path.exists(output_dir):
                os.mkdir(output_dir)
            for amplicon in sample_coverages[sample_name]:
                coverage = sample_coverages[sample_name][amplicon]
                amplicon_file = os.path.join(genome, amplicon + '.fasta')
                read_file = os.path.join(output_dir, amplicon)
                # subprocess_command = 'art_bin_MountRainier/art_illumina --quiet -sam -i ' + amplicon_file + \
                    #    ' --paired -l 250 -f ' + str(coverage) + ' -m 2500 -s 50 -o ' + read_file
                subprocess_command = 'art_bin_MountRainier/art_illumina --quiet -amp -p -sam -na -i ' + amplicon_file + \
                        ' -l 250 -f ' + str(coverage) + ' -o ' + read_file
                subprocess.run(subprocess_command, shell=True, check=True)

        def simulate_badreads(genome,
                              output,
                              sample_coverages):
            sample_name = os.path.basename(genome[:-1])
            output_dir = os.path.join(output, sample_name)
            if not os.path.exists(output_dir):
                os.mkdir(output_dir)
            for amplicon in sample_coverages[sample_name]:
                coverage = sample_coverages[sample_name][amplicon]
                amplicon_file = os.path.join(genome, amplicon + '.fasta')
                read_file = os.path.join(output_dir, amplicon) + '.fastq.gz'
                subprocess_command = 'badread simulate --reference ' + amplicon_file + ' --quantity ' + str(coverage) + 'x \
                |         gzip > ' + read_file
                subprocess.run(subprocess_command, shell=True, check=True)
            return

        # make output dirs
        art_output = os.path.join(output[0], 'ART_output')
        badread_output = os.path.join(output[0], 'Badread_output')
        if not os.path.exists(output[0]):
            os.mkdir(output[0])
        if not os.path.exists(art_output):
            os.mkdir(art_output)
        if not os.path.exists(badread_output):
            os.mkdir(badread_output)
        # list samples and import coverages
        samples = glob.glob(input[0] + '/*/')
        with open(os.path.join(input[0], 'amplicon_coverages.pickle'), 'rb') as coverageHandle:
            sample_coverages = pickle.load(coverageHandle)
        # if specified, assign samples as illumina or nanopore reads
        random.shuffle(samples)
        num_illumina = round(len(samples)*float(params.prop_illumina))
        illumina_samples = samples[:num_illumina]
        nanopore_samples = samples[num_illumina:]
        # parallelise ART
        illumina_list = [
            illumina_samples[i:i + params.processes] for i in range(0, len(illumina_samples), params.processes)
        ]
        for illumina_subset in tqdm(illumina_list):
            Parallel(n_jobs=params.processes, prefer="threads")(delayed(simulate_ART_reads)(genome,
                                                                                            art_output,
                                                                                            sample_coverages) for genome in illumina_subset)
        # parallelise Badread
        nanopore_list = [
            nanopore_samples[i:i + params.processes] for i in range(0, len(nanopore_samples), params.processes)
        ]
        for nanopore_subset in tqdm(nanopore_list):
            Parallel(n_jobs=params.processes, prefer="threads")(delayed(simulate_badreads)(genome,
                                                                                           badread_output,
                                                                                           sample_coverages) for genome in nanopore_subset)

rule cat_reads:
    input:
        rules.simulate_reads.output
    output:
        directory('concatenated_reads')
    run:
        samples = glob.glob(input[0] + '/ART_output/*/')
        if not os.path.exists(output[0]):
            os.mkdir(output[0])
        for genome in samples:
            # concat forward reads
            forward_reads = sorted(glob.glob(os.path.join(genome, '*1.fq')))
            forward_filename = os.path.join(output[0], os.path.basename(genome)) + '_1.fq'
            fw_concat_command = 'cat ' + ' '.join(forward_reads) + ' > ' + forward_filename
            shell(fw_concat_command)
            # concat reverse reads
            reverse_reads = sorted(glob.glob(os.path.join(genome, '*2.fq')))
            reverse_filename = os.path.join(output[0], os.path.basename(genome)) + '_2.fq'
            rv_concat_command = 'cat ' + ' '.join(reverse_reads) + ' > ' + reverse_filename
            shell(rv_concat_command)

rule run_viridian:
    input:
        simulated_genomes=rules.cat_reads.output,
        primer_bed=config['viridian_workflow']['primer_bed'],
        reference_genome=config["reference_genome"]
    output:
        directory('viridian_output')
    params:
        processes=config['processes']
    run:
        def call_vidian_workflow(primer_bed,
                                 reference_genome,
                                 sample,
                                 output):
            fw_read = sample
            rv_read = sample.replace('_1.fq', '_2.fq')
            output_dir = os.path.join(output, os.path.basename(sample).replace('_1.fq', ''))
            viridian_command = 'singularity run viridian_workflow/viridian_workflow.img run_one_sample '
            viridian_command += reference_genome + ' ' + primer_bed + ' ' + fw_read + ' ' + rv_read + ' ' + output_dir + '/'
            subprocess.run(viridian_command, shell=True, check=True)
        samples = glob.glob(os.path.join(input[0], '*_1.fq'))
        if not os.path.exists(output[0]):
            os.mkdir(output[0])
        # parallelise viridian
        job_list = [
            samples[i:i + params.processes] for i in range(0, len(samples), params.processes)
        ]
        for subset in job_list:
            Parallel(n_jobs=params.processes, prefer="threads")(delayed(call_vidian_workflow)(input[1],
                                                                                              input[2],
                                                                                              sample,
                                                                                              output[0]) for sample in subset)

rule assess_assemblies:
    input:
        viridian_assemblies=rules.run_viridian.output,
        simulated_genomes=rules.split_sequences.output
    output:
        directory('viridian_performance')
    params:
        primer_positions=config['split_amplicons']['primer_positions']
    run:
        if not os.path.exists(output[0]):
            os.mkdir(output[0])
        # import bed file as list
        with open(params.primer_positions, "r") as pos:
            primer_metadata = pos.read().splitlines()
        # list viridian-generated assemblies
        assemblies = glob.glob(input[0] + '/*/')
        # get truth set of amplicons from simulated sequences
        mismatch_counts = []
        for assem in assemblies:
            simulated_genome = os.path.join(input[1], os.path.basename(assem[:-1]) + '.fasta')
            truth =  trim_genome(primer_metadata,
                                simulated_genome)
            mismatches = count_errors(truth,
                                    assem)
            mismatch_counts.append(mismatches)
        # visualise distribution of mismatch counts
        counted =  Counter(mismatch_counts)
        pp = PdfPages(os.path.join(output[0], 'viridian_mismatch_counts.pdf'))
        plt.figure()
        plt.hist(x=counted, bins='auto', color='#0504aa', alpha=0.7, rwidth=0.85)
        plt.grid(axis='y', alpha=0.75)
        plt.xlabel('Mismatches')
        plt.ylabel('Frequency')
        plt.title('Distribution of base mismatches in truth genomes vs viridian assemblies')
        pp.savefig()
        pp.close()

rule clean_outputs:
    input:
        VGsim_output=rules.VGsim_tree.output,
        phastSim_output=rules.phastSim_evolution.output,
        split_sequences=rules.split_sequences.output,
        split_amplicons=rules.split_amplicons.output,
        read_output=rules.simulate_reads.output,
        cat_reads=rules.cat_reads.output,
        run_viridian=rules.run_viridian.output,
        viridian_eval=rules.assess_assemblies.output
    shell:
        'rm -rf {input.VGsim_output} {input.phastSim_output} {input.split_sequences} {input.split_amplicons} \
                {input.read_output} {input.cat_reads} {input.run_viridian} {input.viridian_eval}'
