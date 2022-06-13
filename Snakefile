import glob
import gzip
from joblib import Parallel, delayed
import os
import pickle
import pysam
import random
import shutil
from tqdm import tqdm

from scripts.error_modes import clean_genome, find_primer_scheme
from scripts.phylogenies import combine_vcfs, initiate_phylogeny, add_samples, optimise_phylogeny
from scripts.make_plots import run_cte, run_varifier, generate_plots, generate_heatmap, plot_varifier_calls, pairwise_compare, count_dropped_bases

configfile: 'config.yml'

rule all:
    input:
        "phylogeny_comparisons",
        "varifier_statistics",
        "cte_artic_output",
        "cte_viridian_output"

rule VGsim_tree:
    input:
        rate_file=config['VGsim']['rate_file'],
        pp_population_model_file=config['VGsim']['pp_population_model_file'],
        mg_population_model_file=config['VGsim']['mg_population_model_file']
    output:
        'newick_output.nwk'
    singularity:
        'viridian_simulations_workflow.img'
    params:
        seed=config['seed'],
        iterations=config['VGsim']['iterations'],
        sample_size=config['VGsim']['sample_size'],
        container_dir=config["container_directory"]
    shell:
        'singularity run {params.container_dir}/images/VGsim.img {input.rate_file} -it {params.iterations} -pm {input.pp_population_model_file} \
                {input.mg_population_model_file} -seed {params.seed} -nwk --sampleSize {params.sample_size}'

rule phastSim_evolution:
    input:
        tree_file=rules.VGsim_tree.output,
        reference_genome="omicron.fasta"#config["reference_genome"]
    output:
        directory(config['phastSim']['output_dir'])
    params:
        seed=config['seed'],
        container_dir=config["container_directory"],
        substitution_rate=config['phastSim']["substitution_rate"]
    shell:
        'mkdir {output} && singularity run {params.container_dir}/images/phastSim.img --outpath {output}/ --seed {params.seed} --createFasta \
                --createInfo --createNewick --createPhylip --treeFile {input.tree_file} \
                --invariable 0.1 --alpha {params.substitution_rate} --omegaAlpha {params.substitution_rate} --hyperMutProbs 0.001 0.001 --hyperMutRates 2.0 5.0 --codon \
                --reference {input.reference_genome} --createMAT'

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
            with open(os.path.join(output[0], sequence.splitlines()[0] + '_a.fasta'), 'w') as outSeq:
                outSeq.write('>' + sequence)

rule split_amplicons:
    input:
        split_sequences=rules.split_sequences.output,
        phastSim_dir=rules.phastSim_evolution.output
    output:
        directory('amplicon_sequences')
    params:
        primer_scheme=config['primer_scheme'],
        scheme_dir=config["primer_scheme_dir"],
        seed=config["seed"],
        random_dropout_probability=config['split_amplicons']['random_dropout_probability'],
        primer_dimer_probability=config['split_amplicons']['primer_dimer_probability'],
        match_coverage_mean=config['split_amplicons']['match_coverage_mean'],
        match_coverage_sd=config['split_amplicons']['match_coverage_sd'],
        mismatch_coverage_mean=config['split_amplicons']['mismatch_coverage_mean'],
        mismatch_coverage_sd=config['split_amplicons']['mismatch_coverage_sd'],
        container_dir=config["container_directory"]
    threads:
        config['threads']
    shell:
        "python scripts/error_modes.py --scheme {params.primer_scheme} --scheme-dir {params.scheme_dir} --input-dir {input.split_sequences} --output-dir {output} --container-dir {params.container_dir} \
            --seed {params.seed} --dropout-prob {params.random_dropout_probability} --dimer-prob {params.primer_dimer_probability} \
            --match-mean {params.match_coverage_mean} --match-sd {params.match_coverage_sd} --mismatch-mean {params.mismatch_coverage_mean} \
            --mismatch-sd {params.mismatch_coverage_sd} --threads {threads}"

rule simulate_reads:
    input:
        rules.split_amplicons.output
    output:
        directory('read_output')
    params:
        divide_genomes=config["divide_genomes"],
        prop_illumina=config['proportion_illumina'],
        read_length=config["simulate_reads"]["illumina_read_length"],
        seed=config["seed"],
        container_dir=config["container_directory"]
    threads:
        config['threads']
    run:
        def simulate_ART_reads(genome,
                               output,
                               sample_coverages,
                               read_length,
                               container_dir):
            """Function to run ART on amplicon sequences per simulated genomic sequence"""
            sample_name = os.path.basename(genome[:-1])
            output_dir = os.path.join(output, sample_name)
            if not os.path.exists(output_dir):
                os.mkdir(output_dir)
            for amplicon in sample_coverages[sample_name]:
                coverage = sample_coverages[sample_name][amplicon]
                amplicon_file = os.path.join(genome, amplicon + '.fasta')
                read_file = os.path.join(output_dir, amplicon)
                shell_command = 'singularity run ' + container_dir +'/images/ART.img --quiet -amp -p -sam -na -i ' + amplicon_file + \
                        ' -l ' + read_length + ' -f ' + str(coverage) + ' -o ' + read_file
                shell(shell_command)

        def simulate_badreads(genome,
                              output,
                              sample_coverages,
                              container_dir):
            sample_name = os.path.basename(genome[:-1])
            output_dir = os.path.join(output, sample_name)
            if not os.path.exists(output_dir):
                os.mkdir(output_dir)
            for amplicon in sample_coverages[sample_name]:
                coverage = sample_coverages[sample_name][amplicon]
                amplicon_file = os.path.join(genome, amplicon + '.fasta')
                read_file = os.path.join(output_dir, amplicon) + '.fastq.gz'
                # skip badread if the coverage is 0
                if str(coverage) == "0":
                    continue
                shell_command = 'singularity run ' + container_dir + '/images/Badread.img simulate --reference ' + amplicon_file + ' --quantity ' + str(coverage) + 'x \
                |         gzip > ' + read_file
                shell(shell_command)
            return

        # make output dirs
        art_output = os.path.join(output[0], 'ART_output')
        badread_output = os.path.join(output[0], 'Badread_output')
        directories = [output[0], art_output, badread_output]
        for subdir in directories:
            if not os.path.exists(subdir):
                os.mkdir(subdir)
        # list samples and import coverages
        samples = glob.glob(input[0] + '/*/')
        with open(os.path.join(input[0], 'amplicon_coverages.pickle'), 'rb') as coverageHandle:
            sample_coverages = pickle.load(coverageHandle)
        # if specified, assign samples as illumina or nanopore reads
        if params.divide_genomes:
            random.seed(params.seed)
            random.shuffle(samples)
            num_illumina = round(len(samples)*float(params.prop_illumina))
            illumina_samples = samples[:num_illumina]
            nanopore_samples = samples[num_illumina:]
        else:
            illumina_samples = samples
            nanopore_samples = samples
        # parallelise ART
        illumina_list = [
            illumina_samples[i:i + threads] for i in range(0, len(illumina_samples), threads)
        ]
        for illumina_subset in tqdm(illumina_list):
            Parallel(n_jobs=threads, prefer="threads")(delayed(simulate_ART_reads)(genome,
                                                                                    art_output,
                                                                                    sample_coverages,
                                                                                    str(params.read_length),
                                                                                    params.container_dir) for genome in illumina_subset)
        # parallelise Badread
        nanopore_list = [
            nanopore_samples[i:i + threads] for i in range(0, len(nanopore_samples), threads)
        ]
        for nanopore_subset in tqdm(nanopore_list):
            Parallel(n_jobs=threads, prefer="threads")(delayed(simulate_badreads)(genome,
                                                                                    badread_output,
                                                                                    sample_coverages,
                                                                                    params.container_dir) for genome in nanopore_subset)

rule cat_reads:
    input:
        rules.simulate_reads.output
    output:
        directory('concatenated_reads')
    run:
        art_samples = glob.glob(input[0] + '/ART_output/*/')
        badread_samples = glob.glob(input[0] + '/Badread_output/*/')
        if not os.path.exists(output[0]):
            os.mkdir(output[0])
            os.mkdir(os.path.join(output[0], "ART_output"))
            os.mkdir(os.path.join(output[0], "Badread_output"))
        # cat ART reads
        for genome in art_samples:
            # concat forward reads
            forward_reads = sorted(glob.glob(os.path.join(genome, '*1.fq')))
            forward_filename = os.path.join(output[0], 'ART_output', os.path.basename(genome[:-1])) + '_1.fastq'
            fw_concat_command = 'cat ' + ' '.join(forward_reads) + ' > ' + forward_filename
            shell(fw_concat_command)
            # concat reverse reads
            reverse_reads = sorted(glob.glob(os.path.join(genome, '*2.fq')))
            reverse_filename = os.path.join(output[0], 'ART_output', os.path.basename(genome[:-1])) + '_2.fastq'
            rv_concat_command = 'cat ' + ' '.join(reverse_reads) + ' > ' + reverse_filename
            shell(rv_concat_command)
        # cat Badreads
        for sample in badread_samples:
            # concat reads
            reads = sorted(glob.glob(os.path.join(sample, '*.fastq.gz')))
            new_filenames = []
            # gunzip read files for concatenation
            for read_file in reads:
                new_name = read_file.replace(".fastq.gz", ".fastq")
                with gzip.open(read_file, 'rb') as f_in:
                    with open(new_name, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                new_filenames.append(new_name)
            filename = os.path.join(output[0], 'Badread_output', os.path.basename(sample[:-1])) + '.fastq'
            concat_command = 'cat ' + ' '.join(new_filenames) + ' > ' + filename
            shell(concat_command)

rule run_viridian:
    input:
        simulated_reads=rules.cat_reads.output,
        reference_genome=config["reference_genome"]
    output:
        directory('viridian_assemblies')
    params:
        container_dir=config["container_directory"]
    threads:
        config['threads']
    run:
        def illumina_viridian_workflow(reference_genome,
                                       sample,
                                       output,
                                       container_dir):
            """Function to run viridian on ART read sets"""
            fw_read = sample
            rv_read = sample.replace('_1', '_2')
            output_dir = os.path.join(output, os.path.basename(sample).replace('_1.fastq', ''))
            viridian_command = "singularity run " + container_dir + "/viridian_workflow/viridian_workflow.img run_one_sample \
                    --tech illumina \
                    --ref_fasta " + reference_genome + " \
                    --reads1 " + fw_read +" \
                    --keep_bam \
                    --reads2 " + rv_read + " \
                    --outdir " + output_dir + "/"
            shell(viridian_command)

        def nanopore_viridian_workflow(reference_genome,
                                       sample,
                                       output,
                                       container_dir):
            """Function to run viridian on Badread read sets"""
            output_dir = os.path.join(output, os.path.basename(sample).replace('.fastq', ''))
            viridian_command = "singularity run " + container_dir + "/viridian_workflow/viridian_workflow.img run_one_sample \
                    --tech ont \
                    --ref_fasta " + reference_genome + " \
                    --reads " + sample + " \
                    --keep_bam \
                    --outdir " + output_dir + "/"
            shell(viridian_command)

        art_samples = glob.glob(os.path.join(input[0], "ART_output", '*_1.fastq'))
        badread_samples = glob.glob(os.path.join(input[0], "Badread_output", '*.fastq'))
        # make necessary directories
        directories = [output[0], os.path.join(output[0], "ART_assemblies"), os.path.join(output[0], "Badread_assemblies")]
        if not os.path.exists(output[0]):
            for directory in directories:
                os.mkdir(directory)
        # parallelise assembly of ART reads
        art_list = [
            art_samples[i:i + threads] for i in range(0, len(art_samples), threads)
        ]
        for art_set in art_list:
            Parallel(n_jobs=threads, prefer="threads")(delayed(illumina_viridian_workflow)(input[1],
                                                                                           sample,
                                                                                           os.path.join(output[0], "ART_assemblies"),
                                                                                           params.container_dir) for sample in art_set)
        # parallelise assembly of Badread reads
        badread_list = [
            badread_samples[i:i + threads] for i in range(0, len(badread_samples), threads)
        ]
        for bad_set in badread_list:
            Parallel(n_jobs=threads, prefer="threads")(delayed(nanopore_viridian_workflow)(input[1],
                                                                                           sample,
                                                                                           os.path.join(output[0], "Badread_assemblies"),
                                                                                           params.container_dir) for sample in bad_set)

rule mask_assemblies:
    input:
        simulated_genomes=rules.split_sequences.output,
        amplicon_sequences=rules.split_amplicons.output
    output:
        directory("masked_truth_assemblies")
    params:
        mask_assemblies=config['mask_assemblies']["apply_mask"]
    run:
        # make directory
        if not os.path.exists(output[0]):
            os.mkdir(output[0])
        # mask bases outside of the amplicon scheme
        with open(os.path.join(input[1], 'amplicon_statistics.pickle'), 'rb') as statIn:
            amplicon_stats = pickle.load(statIn)
        # mask bases outside of scheme per sample
        for sample in amplicon_stats:
            sample_stats = amplicon_stats[sample]
            amplicons = list(sample_stats.keys())
            first_amplicon_start = amplicon_stats[sample][amplicons[0]]["amplicon_start"]
            last_amplicon_end = amplicon_stats[sample][amplicons[-1]]["amplicon_end"]
            # import truth genome and cut off bases at the ends
            filename = str(sample) + ".fasta"
            sample_name, sample_sequence = clean_genome(os.path.join(input[0], filename))
            sample_sequence = list(sample_sequence)
            sample_sequence = sample_sequence[first_amplicon_start: last_amplicon_end]
            # write out the masked simulated sequence
            with open(os.path.join(output[0], filename), "w") as outGen:
                outGen.write("\n".join([">" + sample_name, "".join(sample_sequence)]))
        # if mask sequences off then just copy the unmasked sequences
        if not params.mask_assemblies:
            pass
        else:
            # iterate through samples and check if amplicons have low coverage, if so then mask their sequence
            for sample in amplicon_stats:
                sample_stats = amplicon_stats[sample]
                all_amplicons = list(sample_stats.keys())
                # import the simulated squence
                filename = str(sample) + ".fasta"
                sample_name, sample_sequence = clean_genome(os.path.join(output[0], filename))
                sample_sequence = list(sample_sequence)
                # look through the amplicon statistics to see if an amplicon needs to be masked
                for amplicon in range(len(all_amplicons)-1):
                    if "primer_SNP" in sample_stats[all_amplicons[amplicon]]["errors"] \
                        or "random_dropout" in sample_stats[all_amplicons[amplicon]]["errors"] \
                        or "primer_dimer" in sample_stats[all_amplicons[amplicon]]["errors"]:
                        # we are masking regions with low coverage that are not covered by the adjacent amplicons
                        if not amplicon == 0:
                            mask_start = sample_stats[all_amplicons[amplicon-1]]["amplicon_end"]
                            #mask_start = sample_stats[all_amplicons[amplicon]]["amplicon_start"]
                        else:
                            mask_start = sample_stats[all_amplicons[amplicon]]["amplicon_start"]
                        mask_end = sample_stats[all_amplicons[amplicon+1]]["amplicon_start"]
                        #mask_end = sample_stats[all_amplicons[amplicon]]["amplicon_end"]
                        # replace the masked sequence with Ns
                        sample_sequence[mask_start:mask_end] = list("N"*(mask_end-mask_start))
                # write out the masked simulated sequence
                with open(os.path.join(output[0], filename), "w") as outGen:
                    outGen.write("\n".join([">" + sample_name, "".join(sample_sequence)]))

rule truth_vcfs:
    input:
        masked_assemblies=rules.mask_assemblies.output,
        reference_genome=config["reference_genome"],
        amplicon_sequences=rules.split_amplicons.output
    output:
        directory("truth_vcfs")
    params:
        container_dir=config["container_directory"],
        primer_scheme=config['primer_scheme']
    threads:
        config['threads']
    run:
        def run_varifier(assembly,
                         covered,
                         dropped_amplicons,
                         primer_df,
                         reference,
                         output_dir,
                         container_dir):
            """Run varifier make_truth_vcf on the masked assemblies"""
            covered_start = covered["start"]
            covered_end = covered["end"]
            varifier_command = "singularity run " + container_dir + "/varifier/varifier.img make_truth_vcf --global_align "
            varifier_command += "--global_align_min_coord " + covered_start + " --global_align_max_coord " + covered_end
            varifier_command += " " + assembly + " " + reference + " " + output_dir
            shell(varifier_command)
            # append dropped amplicon information to the truth vcf
            with open(os.path.join(output_dir, "04.truth.vcf"), "r") as truth_vcf_in:
                truth_vcf = truth_vcf_in.read()
            to_add = []
            variant_count = int(truth_vcf.splitlines()[-1].split("\t")[2])
            for dropped in dropped_amplicons:
                start_pos = primer_df.loc[primer_df['name'] == dropped[0]].reset_index(drop=True)["ref_start"][0]
                end_pos = primer_df.loc[primer_df['name'] == dropped[1]].reset_index(drop=True)["ref_end"][0]
                variant_count += 1
                variant_line = "MN908947.3\t" + str(start_pos) + "\t" + str(variant_count) + "\tG\tN\t.\tDROPPED_AMP\tAMP_START=" + str(int(start_pos)-1) + ";AMP_END=" + str(int(end_pos)-1) + "\tGT\t1/1\n"
                truth_vcf += variant_line
            truth_vcf = truth_vcf.split("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample")
            truth_vcf_variants = truth_vcf[1].splitlines()[1:]
            first_line = [truth_vcf[1].splitlines()[1]]
            truth_vcf_variants.sort(key=lambda x: int(x.split("\t")[1]))
            truth_vcf[1] = "\n".join(first_line + truth_vcf_variants)
            truth_vcf = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample".join(truth_vcf)
            with open(os.path.join(output_dir, "04.truth_dropped.vcf"), "w") as truth_vcf_out:
                truth_vcf_out.write("\n".join(truth_vcf.splitlines()[1:]))

        # make directory
        if not os.path.exists(output[0]):
            os.mkdir(output[0])
        # load list of assemblies
        assemblies = glob.glob(os.path.join(input[0], "*.fasta"))
        # import the amplicon statistics file to extract what parts of the assembly are covered by amplicons
        with open(os.path.join(input[2], 'amplicon_statistics.pickle'), 'rb') as statIn:
            amplicon_stats = pickle.load(statIn)
        regions_covered = {}
        amplicons_dropped = {}
        for sample in amplicon_stats:
            amplicons = list(amplicon_stats[sample].keys())
            regions_covered[sample] = {"start": str(amplicon_stats[sample][amplicons[0]]["amplicon_start"]),
                                       "end": str(amplicon_stats[sample][amplicons[len(amplicons)-1]]["amplicon_end"])}
            # record amplicons dropped for each sample so dropped amplicons can be marked in the truth vcf
            amplicons_dropped[sample] = []
            for a in amplicons:
                if amplicon_stats[sample][a]["has_error"] and any(amplicon_stats[sample][a]["errors"][0] == mode for mode in ["primer_SNP", "random_dropout", "primer_dimer"]):
                    amplicons_dropped[sample].append((a.split("---")[0], a.split("---")[1]))
        # import the primer scheme df to get amplicon positions relative to ref
        primer_df, pool1_primers, pool2_primers = find_primer_scheme(params.primer_scheme,
                                                                    "primer_schemes")
        # parallelise make_truth_vcf
        subsetted_assemblies = [
            assemblies[i:i + threads] for i in range(0, len(assemblies), threads)
        ]
        for subset in subsetted_assemblies:
            Parallel(n_jobs=threads, prefer="threads")(delayed(run_varifier)(sample,
                                                                            regions_covered[os.path.basename(sample).replace(".fasta", "")],
                                                                            amplicons_dropped[os.path.basename(sample).replace(".fasta", "")],
                                                                            primer_df,
                                                                            input[1],
                                                                            os.path.join(output[0], os.path.basename(sample).replace(".fasta", "")),
                                                                            params.container_dir) for sample in subset)

rule artic_assemble:
    input:
        simulated_reads=rules.cat_reads.output,
        amplicon_sequences=rules.split_amplicons.output
    output:
        directory("artic_assemblies")
    threads:
        config["threads"]
    params:
        scheme_dir=config["primer_scheme_dir"],
        primer_scheme=config['primer_scheme'],
        illumina_container=config["artic_assemble"]["illumina_workflow_container"],
        nanopore_container=config["artic_assemble"]["nanopore_workflow_container"],
        container_dir=config["container_directory"],
        nextflow_path=config["artic_assemble"]["nextflow_path"]
    run:

        def illumina_artic_assemble(forward_rd,
                                    reverse_rd,
                                    sif_file,
                                    main_nf,
                                    scheme_url,
                                    output_dir,
                                    nextflow_path,
                                    primer_scheme,
                                    regions_covered):
            """run illumina artic nextflow pipeline"""
            if not ".gz" in forward_rd:
                shell("gzip " + forward_rd + " " + reverse_rd)
                forward_rd = forward_rd + ".gz"
                reverse_rd = reverse_rd + ".gz"
            else:
                forward_rd = forward_rd
                reverse_rd = reverse_rd
            shell_command = "venv/bin/python3.9 scripts/run_connor_pipeline.py --sif " + sif_file + " "
            shell_command += "--main_nf " + main_nf + " --outdir " + output_dir.replace(".gz", "") + " "
            shell_command += "--scheme_url " + scheme_url + " --scheme_version " + primer_scheme + " "
            shell_command += "--ilm1 " + forward_rd + " --ilm2 " + reverse_rd + " --nextflow_path " + nextflow_path + " "
            shell_command += "--sample_name " + os.path.basename(forward_rd).replace("_1.fastq", "").replace(".gz", "") + " && "
            if not ".gz" in forward_rd:
                forward_rd = forward_rd + ".gz"
                reverse_rd = reverse_rd + ".gz"
            shell_command += "gunzip " + forward_rd + " " + reverse_rd
            shell(shell_command)
            # cut off ends of assemblies that lie outside of the amplicon scheme
            filename = os.path.join(output_dir.replace(".gz", ""), "consensus.fa")
            sample_name, sample_sequence = clean_genome(filename)
            sample_sequence = list(sample_sequence)
            sample_sequence = sample_sequence[regions_covered["scheme_start"]: regions_covered["scheme_end"]]
            # write out the masked simulated sequence
            with open(filename.replace("consensus", "consensus_trimmed"), "w") as outGen:
                outGen.write("\n".join([">" + sample_name, "".join(sample_sequence)]))

        def nanopore_artic_assemble(read_file,
                                    sif_file,
                                    main_nf,
                                    scheme_url,
                                    output_dir,
                                    nextflow_path,
                                    primer_scheme,
                                    regions_covered):
            """run nanopore artic nextflow pipeline"""
            if not ".gz" in read_file:
                shell("gzip " + read_file)
                read_file = read_file + ".gz"
            else:
                read_file = read_file
            shell_command = "venv/bin/python3.9 scripts/run_connor_pipeline.py --sif " + sif_file + " "
            shell_command += "--main_nf " + main_nf + " --outdir " + output_dir.replace(".gz", "") + " "
            shell_command += "--ont " + read_file + " "
            shell_command += "--scheme_url " + scheme_url + " --scheme_version " + primer_scheme + " --nextflow_path " + nextflow_path + " "
            shell_command += "--sample_name " + os.path.basename(read_file).replace(".fastq", "").replace(".gz", "") + " && "
            if not ".gz" in read_file:
                read_file = read_file + ".gz"
            shell_command += "gunzip " + read_file
            shell(shell_command)
            # cut off ends of assemblies that lie outside of the amplicon scheme
            filename = os.path.join(output_dir.replace(".gz", ""), "consensus.fa")
            sample_name, sample_sequence = clean_genome(filename)
            sample_sequence = list(sample_sequence)
            sample_sequence = sample_sequence[regions_covered["scheme_start"]: regions_covered["scheme_end"]]
            # write out the masked simulated sequence
            with open(filename.replace("consensus", "consensus_trimmed"), "w") as outGen:
                outGen.write("\n".join([">" + sample_name, "".join(sample_sequence)]))

        # list read sets
        art_samples = glob.glob(os.path.join(input[0], "ART_output", '*_1.fastq'))
        art_samples += glob.glob(os.path.join(input[0], "ART_output", '*_1.fastq.gz'))
        badread_samples = glob.glob(os.path.join(input[0], "Badread_output", '*.fastq'))
        badread_samples += glob.glob(os.path.join(input[0], "Badread_output", '*.fastq.gz'))
        # make output directories
        for subdir in [output[0], os.path.join(output[0], "ART_assemblies"), os.path.join(output[0], "Badread_assemblies")]:
            if not os.path.exists(subdir):
                os.mkdir(subdir)
        # run artic pipeline on read sets
        job_threads = int(-(-threads // 2))
        subsetted_art = [
            art_samples[i:i + job_threads] for i in range(0, len(art_samples), job_threads)
        ]
        subsetted_badread = [
            badread_samples[i:i + job_threads] for i in range(0, len(badread_samples), job_threads)
        ]
        # import the primer scheme df to get amplicon scheme start and end relative to ref
        primer_df, pool1_primers, pool2_primers = find_primer_scheme(params.primer_scheme,
                                                                    "primer_schemes")
        with open(os.path.join(input[1], 'amplicon_statistics.pickle'), 'rb') as statIn:
            amplicon_stats = pickle.load(statIn)
        regions_covered = {}
        for sample in amplicon_stats:
            amplicons = list(amplicon_stats[sample].keys())
            scheme_start = int(primer_df.loc[primer_df['name'] == amplicons[0].split("---")[0]].reset_index(drop=True)["ref_start"][0])
            scheme_end = int(primer_df.loc[primer_df['name'] == amplicons[len(amplicons)-1].split("---")[1]].reset_index(drop=True)["ref_start"][0]) + 1
            regions_covered[sample] = {"scheme_start": scheme_start,
                                       "scheme_end": scheme_end}
    
        for subset in tqdm(subsetted_badread):
            Parallel(n_jobs=job_threads, prefer="threads")(delayed(nanopore_artic_assemble)(read,
                                                                                            params.nanopore_container,
                                                                                            os.path.join(params.container_dir, "ncov2019-artic-nf", "main.nf"),
                                                                                            params.scheme_dir,
                                                                                            os.path.join(output[0], "Badread_assemblies", os.path.basename(read).replace(".fastq", "")),
                                                                                            params.nextflow_path,
                                                                                            params.primer_scheme,
                                                                                            regions_covered[os.path.basename(read).replace(".fastq", "").replace(".gz", "")]) for read in subset)
        for subset in tqdm(subsetted_art):
            Parallel(n_jobs=job_threads, prefer="threads")(delayed(illumina_artic_assemble)(read,
                                                                                            read.replace("_1.fastq", "_2.fastq"),
                                                                                            params.illumina_container,
                                                                                            os.path.join(params.container_dir, "ncov2019-artic-nf", "main.nf"),
                                                                                            params.scheme_dir,
                                                                                            os.path.join(output[0], "ART_assemblies", os.path.basename(read).replace("_1.fastq", "")),
                                                                                            params.nextflow_path,
                                                                                            params.primer_scheme,
                                                                                            regions_covered[os.path.basename(read).replace("_1.fastq", "").replace(".gz", "")]) for read in subset)

rule viridian_covid_truth_eval:
    input:
        reference_genome=config["reference_genome"],
        viridian_assemblies=rules.run_viridian.output,
        truth_vcf=rules.truth_vcfs.output,
    output:
        directory("cte_viridian_output")
    threads:
        config["threads"]
    params:
        container_dir=config["container_directory"],
        primer_scheme=config["primer_scheme"]
    run:
        # make output dirs
        for sub_dir in [output[0], os.path.join(output[0], "ART_assemblies"), os.path.join(output[0], "Badread_assemblies")]:
            if not os.path.exists(sub_dir):
                os.mkdir(sub_dir)
        # list viridian assemblies
        art_assemblies = glob.glob(os.path.join(input[1], "ART_assemblies", "*"))
        badread_assemblies = glob.glob(os.path.join(input[1], "Badread_assemblies", "*"))
        # run covid truth eval
        run_cte(params.primer_scheme,
                art_assemblies,
                os.path.join(output[0], "ART_assemblies"),
                input[2],
                params.container_dir,
                "viridian")
        run_cte(params.primer_scheme,
                badread_assemblies,
                os.path.join(output[0], "Badread_assemblies"),
                input[2],
                params.container_dir,
                "viridian")

rule artic_covid_truth_eval:
    input:
        reference_genome=config["reference_genome"],
        artic_assemblies=rules.artic_assemble.output,
        truth_vcf=rules.truth_vcfs.output,
    output:
        directory("cte_artic_output")
    threads:
        config["threads"]
    params:
        container_dir=config["container_directory"],
        primer_scheme=config["primer_scheme"]
    run:
       # make output dirs
        for sub_dir in [output[0], os.path.join(output[0], "ART_assemblies"), os.path.join(output[0], "Badread_assemblies")]:
            if not os.path.exists(sub_dir):
                os.mkdir(sub_dir)
        # list viridian assemblies
        art_assemblies = glob.glob(os.path.join(input[1], "ART_assemblies", "*"))
        badread_assemblies = glob.glob(os.path.join(input[1], "Badread_assemblies", "*"))
        # run covid truth eval
        run_cte(params.primer_scheme,
                art_assemblies,
                os.path.join(output[0], "ART_assemblies"),
                input[2],
                params.container_dir,
                "artic")
        run_cte(params.primer_scheme,
                badread_assemblies,
                os.path.join(output[0], "Badread_assemblies"),
                input[2],
                params.container_dir,
                "artic")

rule make_artic_vcfs:
    input:
        artic_assemblies=rules.artic_assemble.output,
        reference_genome=config["reference_genome"],
        amplicon_sequences=rules.split_amplicons.output
    output:
        directory("artic_vcfs")
    params:
        container_dir=config["container_directory"]
    threads:
        config['threads']
    run:
        def run_varifier(assembly,
                         covered,
                         reference,
                         output_dir,
                         container_dir):
            """Run varifier make_truth_vcf on the masked assemblies"""
            covered_start = covered["start"]
            covered_end = covered["end"]
            varifier_command = "singularity run " + container_dir + "/varifier/varifier.img make_truth_vcf --global_align "
            varifier_command += "--global_align_min_coord " + covered_start + " --global_align_max_coord " + covered_end
            varifier_command += " " + os.path.join(assembly, "consensus_trimmed.fa") + " " + reference + " " + output_dir
            shell(varifier_command)

        # make directory
        for sub_dir in [output[0], os.path.join(output[0], "ART_assemblies"), os.path.join(output[0], "Badread_assemblies")]:
            if not os.path.exists(sub_dir):
                os.mkdir(sub_dir)
        # load list of assemblies
        art_assemblies = glob.glob(os.path.join(input[0], "ART_assemblies", "*"))
        badread_assemblies = glob.glob(os.path.join(input[0], "Badread_assemblies", "*"))
        # import the amplicon statistics file to extract what parts of the assembly are covered by amplicons
        with open(os.path.join(input[2], 'amplicon_statistics.pickle'), 'rb') as statIn:
            amplicon_stats = pickle.load(statIn)
        regions_covered = {}
        for sample in amplicon_stats:
            amplicons = list(amplicon_stats[sample].keys())
            regions_covered[sample] = {"start": str(amplicon_stats[sample][amplicons[0]]["amplicon_start"]),
                                       "end": str(amplicon_stats[sample][amplicons[len(amplicons)-1]]["amplicon_end"])}
        # parallelise make_truth_vcf
        subsetted_art = [
            art_assemblies[i:i + threads] for i in range(0, len(art_assemblies), threads)
        ]
        for subset in subsetted_art:
            Parallel(n_jobs=threads, prefer="threads")(delayed(run_varifier)(sample,
                                                                            regions_covered[os.path.basename(sample)],
                                                                            input[1],
                                                                            os.path.join(output[0], "ART_assemblies", os.path.basename(sample)),
                                                                            params.container_dir) for sample in subset)
        subsetted_badread = [
            badread_assemblies[i:i + threads] for i in range(0, len(badread_assemblies), threads)
        ]
        for subset in subsetted_badread:
            Parallel(n_jobs=threads, prefer="threads")(delayed(run_varifier)(sample,
                                                                            regions_covered[os.path.basename(sample)],
                                                                            input[1],
                                                                            os.path.join(output[0], "Badread_assemblies", os.path.basename(sample)),
                                                                            params.container_dir) for sample in subset)

rule assess_assemblies:
    input:
        viridian_assemblies=rules.run_viridian.output,
        simulated_genomes=rules.mask_assemblies.output,
        reference_genome=config["reference_genome"],
        truth_vcf=rules.truth_vcfs.output,
        split_amplicons=rules.split_amplicons.output,
        artic_vcfs=rules.make_artic_vcfs.output,
        artic_assemblies=rules.artic_assemble.output,
        unmasked_genomes=rules.split_sequences.output,
        artic_cte=rules.artic_covid_truth_eval.output,
        viridian_cte=rules.viridian_covid_truth_eval.output
    params:
        container_dir=config["container_directory"],
        scheme=config['primer_scheme']
    threads:
        config['threads']
    output:
        directory('varifier_statistics')
    run:
        directories = [output[0], os.path.join(output[0], "viridian_assemblies"), os.path.join(output[0], "artic_assemblies"), \
            os.path.join(output[0], "viridian_assemblies", "ART_assemblies"), os.path.join(output[0], "viridian_assemblies", "Badread_assemblies"), \
            os.path.join(output[0], "artic_assemblies", "ART_assemblies"), os.path.join(output[0], "artic_assemblies", "Badread_assemblies")]
        if not os.path.exists(output[0]):
            for folder in directories:
                os.mkdir(folder)
        art_assemblies = glob.glob(input[0] + '/ART_assemblies/*')
        badread_assemblies = glob.glob(input[0] + '/Badread_assemblies/*')
        # import amplicon statistics
        with open(os.path.join(input[4], 'amplicon_statistics.pickle'), 'rb') as statIn:
            amplicon_stats = pickle.load(statIn)
        # iterate through assemblies, run varifier and generate plots 
        viridian_art_snps = generate_plots(art_assemblies,
                                        amplicon_stats,
                                        input[3],
                                        input[2],
                                        input[1],
                                        "ART",
                                        threads,
                                        output[0],
                                        params.container_dir)
        viridian_badread_snps = generate_plots(badread_assemblies,
                                            amplicon_stats,
                                            input[3],
                                            input[2],
                                            input[1],
                                            "Badread",
                                            threads,
                                            output[0],
                                            params.container_dir)
        # generate plots for the artic assemblies
        artic_art_assemblies = glob.glob(input[6] + '/ART_assemblies/*')
        artic_badread_assemblies = glob.glob(input[6] + '/Badread_assemblies/*')
        artic_art_snps = generate_plots(artic_art_assemblies,
                                            amplicon_stats,
                                            input[3],
                                            input[2],
                                            input[1],
                                            "artic_ART",
                                            threads,
                                            output[0],
                                            params.container_dir,
                                            input[5],
                                            os.path.join(input[0], "ART_assemblies"))
        artic_badread_snps = generate_plots(artic_badread_assemblies,
                                            amplicon_stats,
                                            input[3],
                                            input[2],
                                            input[1],
                                            "artic_Badread",
                                            threads,
                                            output[0],
                                            params.container_dir,
                                            input[5],
                                            os.path.join(input[0], "Badread_assemblies"))
        # let's see what happens if we assess SNPs called by the assemblers instead of truth SNPs
        plot_varifier_calls(viridian_art_snps,
                            viridian_badread_snps,
                            artic_art_snps,
                            artic_badread_snps,
                            output[0])
        # generate heatmap plots for the base by base calls by cte
                # generate plot for results of covid truth eval output
        generate_heatmap(os.path.join(input[8], "ART_assemblies", "OUT", "Processing"), "ART", output[0])
        generate_heatmap(os.path.join(input[8], "Badread_assemblies", "OUT", "Processing"), "Badread", output[0])
        generate_heatmap(os.path.join(input[9], "ART_assemblies", "OUT", "Processing"), "ART", output[0])
        generate_heatmap(os.path.join(input[9], "Badread_assemblies", "OUT", "Processing"), "Badread", output[0])
        # generate matrices of pairwise SNP comparisons in viridian vs artic assemblies
        pairwise_compare(art_assemblies,
                         os.path.join(input[6], "ART_assemblies"),
                         output[0],
                         "viridian_artic",
                         params.container_dir)
        pairwise_compare(artic_art_assemblies,
                         input[7],
                         output[0],
                         "artic_simulated",
                         params.container_dir)
        pairwise_compare(art_assemblies,
                         input[7],
                         output[0],
                         "viridian_simulated",
                         params.container_dir)
        # count the actual number of low coverage bases
        count_dropped_bases(input[4],
                            output[0])


rule build_viridian_phylogeny:
    input:
        viridian_assemblies=rules.run_viridian.output,
        reference_genome=config["reference_genome"]
    output:
        directory("viridian_phylogenies")
    threads:
        config["threads"]
    params:
        batch_size=config["phylogenies"]["batch_size"]
    run:
        # list assemblies
        art_assemblies = glob.glob(os.path.join(input[0], "ART_assemblies", "*"))
        badread_assemblies = glob.glob(os.path.join(input[0], "Badread_assemblies", "*"))
        # order assemblies so ART and Badread assemblies are added to the tree in the same order
        sorted(art_assemblies, key=str.lower)
        sorted(badread_assemblies, key=str.lower)
        # concatenate assemblies into single file for alignment
        art_fasta = []
        for assem in art_assemblies:
            with open(os.path.join(assem, "consensus.fa"), "r") as inAssem:
                art_fasta.append(">" + os.path.basename(assem) + "\n" + "".join(inAssem.read().splitlines()[1:]))
        # make output directory if it does not already exist
        for directory in [output[0], os.path.join(output[0], "ART_assemblies"), os.path.join(output[0], "Badread_assemblies")]:
            if not os.path.exists(directory):
                os.mkdir(directory)
        # build ART and Badread trees using USHER
        for method in [art_assemblies, badread_assemblies]:
            if method == art_assemblies:
                out_dir = os.path.join(output[0], "ART_assemblies")
            if method == badread_assemblies:
                out_dir =  os.path.join(output[0], "Badread_assemblies")
            # make vcf files to add in batches
            batched_files = [
                method[i:i + params.batch_size] for i in range(0, len(method), params.batch_size)
            ]
            out_vcfs = []
            for batch in tqdm(batched_files):
                out_vcfs.append(combine_vcfs(batch,
                                out_dir,
                                "viridian",
                                input[0]))
            # make trees with one sample to start the phylogeny
            if method == art_assemblies:
                vcf_dir = os.path.join(input[0], "ART_assemblies")
                tree_dir, tree_file = initiate_phylogeny("ART",
                                                         method,
                                                         input[0],
                                                         vcf_dir,
                                                         "viridian",
                                                         output[0],
                                                         input[1])
            else:
                vcf_dir = os.path.join(input[0], "Badread_assemblies")
                tree_dir, tree_file = initiate_phylogeny("Badread",
                                                         method,
                                                         input[0],
                                                         vcf_dir,
                                                         "viridian",
                                                         output[0],
                                                         input[1])
            for combined in out_vcfs:
                # add the samples to the phylogeny in batches
                add_samples(combined,
                            tree_file,
                            tree_dir,
                            threads)
                # optimise the tree after every batch is added
                optimise_phylogeny(tree_file,
                                   threads)

rule build_simulated_phylogeny:
    input:
        simulated_assemblies=rules.mask_assemblies.output,
        viridian_assemblies=rules.run_viridian.output,
        truth_vcfs=rules.truth_vcfs.output,
        reference_genome=config["reference_genome"]
    output:
        directory("simulated_phylogenies")
    threads:
        config["threads"]
    params:
        batch_size=config["phylogenies"]["batch_size"]
    run:
        # list assemblies
        viridian_art_assemblies = glob.glob(os.path.join(input[1], "ART_assemblies", "*"))
        viridian_badread_assemblies = glob.glob(os.path.join(input[1], "Badread_assemblies", "*"))
        # order assemblies so ART and Badread assemblies are added to the tree in the same order
        sorted(viridian_art_assemblies, key=str.lower)
        sorted(viridian_badread_assemblies, key=str.lower)
        # need to assign the simulated assemblies to either ART or Badread, depending on the viridian assembly assignments
        art_assemblies = []
        for a in viridian_art_assemblies:
            assem = os.path.join(input[0], os.path.basename(a))
            art_assemblies.append(assem)
        badread_assemblies = []
        for a in viridian_badread_assemblies:
            assem = os.path.join(input[0], os.path.basename(a))
            badread_assemblies.append(assem)
        # concatenate assemblies into single file for alignment
        art_fasta = []
        for assem in art_assemblies:
            with open(assem + ".fasta", "r") as inAssem:
                art_fasta.append(">" + os.path.basename(assem) + "\n" + "".join(inAssem.read().splitlines()[1:]))
        # make output directory if it does not already exist
        for directory in [output[0], os.path.join(output[0], "ART_assemblies"), os.path.join(output[0], "Badread_assemblies")]:
            if not os.path.exists(directory):
                os.mkdir(directory)
        # build ART and Badread trees using USHER
        for method in [art_assemblies, badread_assemblies]:
            if method == art_assemblies:
                out_dir = os.path.join(output[0], "ART_assemblies")
            if method == badread_assemblies:
                out_dir =  os.path.join(output[0], "Badread_assemblies")
            # make vcf files to add in batches
            batched_files = [
                method[i:i + params.batch_size] for i in range(0, len(method), params.batch_size)
            ]
            out_vcfs = []
            out_vcfs += Parallel(n_jobs=threads, prefer="threads")(delayed(combine_vcfs)(batch,
                                                                                        out_dir,
                                                                                        "simulated",
                                                                                        input[0],
                                                                                        input[2]) for batch in batched_files)
            # make trees with one sample to start the phylogeny
            vcf_dir = input[2]
            if method == art_assemblies:
                tree_dir, tree_file = initiate_phylogeny("ART",
                                                         method,
                                                         input[0],
                                                         vcf_dir,
                                                         "simulated",
                                                         output[0],
                                                         input[3])
            else:
                tree_dir, tree_file = initiate_phylogeny("Badread",
                                                         method,
                                                         input[0],
                                                         vcf_dir,
                                                         "simulated",
                                                         output[0],
                                                         input[3])
            for combined in out_vcfs:
                # add the samples to the phylogeny in batches
                add_samples(combined,
                            tree_file,
                            tree_dir,
                            threads)
                # optimise the tree after every batch is added
                optimise_phylogeny(tree_file,
                                    threads)

rule build_artic_phylogeny:
    input:
        simulated_assemblies=rules.mask_assemblies.output,
        viridian_assemblies=rules.run_viridian.output,
        artic_vcfs=rules.make_artic_vcfs.output,
        reference_genome=config["reference_genome"],
        artic_assemblies=rules.artic_assemble.output
    output:
        directory("artic_phylogenies")
    threads:
        config["threads"]
    params:
        batch_size=config["phylogenies"]["batch_size"]
    run:
        # list assemblies
        viridian_art_assemblies = glob.glob(os.path.join(input[1], "ART_assemblies", "*"))
        viridian_badread_assemblies = glob.glob(os.path.join(input[1], "Badread_assemblies", "*"))
        # order assemblies so ART and Badread assemblies are added to the tree in the same order
        sorted(viridian_art_assemblies, key=str.lower)
        sorted(viridian_badread_assemblies, key=str.lower)
        # need to assign the simulated assemblies to either ART or Badread, depending on the viridian assembly assignments
        art_assemblies = []
        for a in viridian_art_assemblies:
            assem = os.path.join(input[0], os.path.basename(a))
            art_assemblies.append(assem)
        badread_assemblies = []
        for a in viridian_badread_assemblies:
            assem = os.path.join(input[0], os.path.basename(a))
            badread_assemblies.append(assem)
        # concatenate assemblies into single file for alignment
        art_fasta = []
        for assem in art_assemblies:
            with open(assem + ".fasta", "r") as inAssem:
                art_fasta.append(">" + os.path.basename(assem) + "\n" + "".join(inAssem.read().splitlines()[1:]))
        # make output directory if it does not already exist
        for directory in [output[0], os.path.join(output[0], "ART_assemblies"), os.path.join(output[0], "Badread_assemblies")]:
            if not os.path.exists(directory):
                os.mkdir(directory)
        # build ART and Badread trees using USHER
        for method in ["ART_assemblies", "Badread_assemblies"]:
            if method == "ART_assemblies":
                out_dir = os.path.join(output[0], "ART_assemblies")
                in_files = art_assemblies
            if method == "Badread_assemblies":
                out_dir =  os.path.join(output[0], "Badread_assemblies")
                in_files = badread_assemblies
            # make vcf files to add in batches
            batched_files = [
                in_files[i:i + params.batch_size] for i in range(0, len(in_files), params.batch_size)
            ]
            out_vcfs = []
            out_vcfs += Parallel(n_jobs=threads, prefer="threads")(delayed(combine_vcfs)(batch,
                                                                                        out_dir,
                                                                                        "artic",
                                                                                        input[4],
                                                                                        input[2]) for batch in batched_files)
            # make trees with one sample to start the phylogeny
            vcf_dir = input[2]
            if method == "ART_assemblies":
                tree_dir, tree_file = initiate_phylogeny("ART",
                                                         in_files,
                                                         input[0],
                                                         vcf_dir,
                                                         "artic",
                                                         output[0],
                                                         input[3])
            else:
                tree_dir, tree_file = initiate_phylogeny("Badread",
                                                         in_files,
                                                         input[0],
                                                         vcf_dir,
                                                         "artic",
                                                         output[0],
                                                         input[3])
            for combined in out_vcfs:
                # add the samples to the phylogeny in batches
                add_samples(combined,
                            tree_file,
                            tree_dir,
                            threads)
                # optimise the tree after every batch is added
                optimise_phylogeny(tree_file,
                                    threads)

rule compare_phylogenies:
    input:
        viridian_phylogeny=rules.build_viridian_phylogeny.output,
        simulated_phylogeny=rules.build_simulated_phylogeny.output,
        artic_phylogeny=rules.build_artic_phylogeny.output
    output:
        directory("phylogeny_comparisons")
    threads:
        config["threads"]
    run:
        sources = [input[0], input[1], input[2]]
        methods = ["ART", "Badread"]
        for source in sources:
            for directory in [output[0], os.path.join(output[0], source), os.path.join(output[0], source, "ART_assemblies"), os.path.join(output[0], source, "Badread_assemblies")]:
                if not os.path.exists(directory):
                    os.mkdir(directory)
            for method in methods:
                output_tsv = os.path.join(output[0], source, method + "_assemblies", "output.tsv")
                phylogeny = os.path.join(source, method + "_assemblies", method + "_phylo.pb")
                shell("singularity exec singularity/usher/usher.sif matUtils extract --closest-relatives " + output_tsv + " --break-ties --threads {threads} -i " + phylogeny)