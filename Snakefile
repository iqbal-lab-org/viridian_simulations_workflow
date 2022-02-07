import glob
import gzip
from joblib import Parallel, delayed
import os
import pickle
import pysam
import random
import shutil
import subprocess
from tqdm import tqdm

from scripts.error_modes import clean_genome

configfile: 'config.yml'

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
        container_dir=config["container_directory"]
    shell:
        'mkdir {output} && singularity run {params.container_dir}/images/phastSim.img --outpath {output}/ --seed {params.seed} --createFasta \
                --createInfo --createNewick --createPhylip --treeFile {input.tree_file} \
                --invariable 0.1 --alpha 0.0002 --omegaAlpha 0.0002 --hyperMutProbs 0.001 0.001 --hyperMutRates 2.0 5.0 --codon \
                --reference {input.reference_genome}'

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
        split_sequences=rules.split_sequences.output,
        phastSim_dir=rules.phastSim_evolution.output
    output:
        directory('amplicon_sequences')
    params:
        primer_scheme=config['primer_scheme'],
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
        "python scripts/error_modes.py --scheme {params.primer_scheme} --input-dir {input.split_sequences} --output-dir {output} --container-dir {params.container_dir} \
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
                subprocess_command = 'singularity run ' + container_dir +'/images/ART.img --quiet -amp -p -sam -na -i ' + amplicon_file + \
                        ' -l ' + read_length + ' -f ' + str(coverage) + ' -o ' + read_file
                subprocess.run(subprocess_command, shell=True, check=True)

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
                subprocess_command = 'singularity run ' + container_dir + '/images/Badread.img simulate --reference ' + amplicon_file + ' --quantity ' + str(coverage) + 'x \
                |         gzip > ' + read_file
                subprocess.run(subprocess_command, shell=True, check=True)
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
            forward_filename = os.path.join(output[0], 'ART_output', os.path.basename(genome[:-1])) + '_1.fq'
            fw_concat_command = 'cat ' + ' '.join(forward_reads) + ' > ' + forward_filename
            shell(fw_concat_command)
            # concat reverse reads
            reverse_reads = sorted(glob.glob(os.path.join(genome, '*2.fq')))
            reverse_filename = os.path.join(output[0], 'ART_output', os.path.basename(genome[:-1])) + '_2.fq'
            rv_concat_command = 'cat ' + ' '.join(reverse_reads) + ' > ' + reverse_filename
            shell(rv_concat_command)
        # cat Badreads
        for sample in badread_samples:
            # concat reads
            reads = sorted(glob.glob(os.path.join(sample, '*.fastq.gz')))
            new_filenames = []
            # gunzip read files for concatenation
            for read_file in reads:
                new_name = read_file.replace(".fastq.gz", ".fq")
                with gzip.open(read_file, 'rb') as f_in:
                    with open(new_name, 'wb') as f_out:
                        shutil.copyfileobj(f_in, f_out)
                new_filenames.append(new_name)
            filename = os.path.join(output[0], 'Badread_output', os.path.basename(sample[:-1])) + '.fq'
            concat_command = 'cat ' + ' '.join(new_filenames) + ' > ' + filename
            shell(concat_command)

rule run_viridian:
    input:
        simulated_genomes=rules.cat_reads.output,
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
            output_dir = os.path.join(output, os.path.basename(sample).replace('_1.fq', ''))
            viridian_command = "singularity run " + container_dir + "/viridian_workflow/viridian_workflow.img run_one_sample \
                    --tech illumina \
                    --ref_fasta " + reference_genome + " \
                    --reads1 " + fw_read +" \
                    --keep_bam \
                    --reads2 " + rv_read + " \
                    --outdir " + output_dir + "/"
            subprocess.run(viridian_command, shell=True, check=True)

        def nanopore_viridian_workflow(reference_genome,
                                       sample,
                                       output,
                                       container_dir):
            """Function to run viridian on Badread read sets"""
            output_dir = os.path.join(output, os.path.basename(sample).replace('.fq', ''))
            viridian_command = "singularity run " + container_dir + "/viridian_workflow/viridian_workflow.img run_one_sample \
                    --tech ont \
                    --ref_fasta " + reference_genome + " \
                    --reads " + sample + " \
                    --keep_bam \
                    --outdir " + output_dir + "/"
            subprocess.run(viridian_command, shell=True, check=True)

        art_samples = glob.glob(os.path.join(input[0], "ART_output", '*_1.fq'))
        badread_samples = glob.glob(os.path.join(input[0], "Badread_output", '*.fq'))
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
        # if mask sequences off then just copy the unmasked sequences
        if not params.mask_assemblies:
            shell("cp -a /" + input[0] + "/. /" + output[0] + "/")
        else:
            # import the amplicon statistics file to see which amplicons to mask
            with open(os.path.join(input[1], 'amplicon_statistics.pickle'), 'rb') as statIn:
                amplicon_stats = pickle.load(statIn)
            # iterate through samples and check if amplicons have low coverage, if so then mask their sequence
            for sample in amplicon_stats:
                sample_stats = amplicon_stats[sample]
                all_amplicons = list(sample_stats.keys())[:-1]
                # import the simulated squence
                filename = str(sample) + ".fasta"
                sample_name, sample_sequence = clean_genome(os.path.join(input[0], filename))
                sample_sequence = list(sample_sequence)
                # look through the amplicon statistics to see if an amplicon needs to be masked
                for amplicon in range(len(all_amplicons)-1):
                    if "primer_SNP" in sample_stats[all_amplicons[amplicon]]["errors"] \
                        or "random_dropout" in sample_stats[all_amplicons[amplicon]]["errors"] \
                        or "primer_dimer" in sample_stats[all_amplicons[amplicon]]["errors"]:
                        # we are masking regions with low coverage that are not covered by the adjacent amplicons
                        mask_start = sample_stats[all_amplicons[amplicon-1]]["right_primer_start"]
                        mask_end = sample_stats[all_amplicons[amplicon+1]]["left_primer_end"]
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
            #varifier_command = 'singularity exec "docker://quay.io/iqballab/varifier" varifier make_truth_vcf --global_align '
            varifier_command += "--global_align_min_coord " + covered_start + " --global_align_max_coord " + covered_end
            varifier_command += " " + assembly + " " + reference + " " + output_dir
            subprocess.run(varifier_command, shell=True, check=True)

        # make directory
        if not os.path.exists(output[0]):
            os.mkdir(output[0])
        # load list of assemblies
        assemblies = glob.glob(os.path.join(input[0], "*.fasta"))
        # import the amplicon statistics file to extract what parts of the assembly are covered by amplicons
        with open(os.path.join(input[2], 'amplicon_statistics.pickle'), 'rb') as statIn:
            amplicon_stats = pickle.load(statIn)
        regions_covered = {}
        for sample in amplicon_stats:
            amplicons = list(amplicon_stats[sample].keys())[:-1]
            regions_covered[sample] = {"start": str(amplicon_stats[sample][amplicons[0]]["amplicon_start"]),
                                       "end": str(amplicon_stats[sample][amplicons[len(amplicons)-1]]["amplicon_end"])}
        # parallelise make_truth_vcf
        subsetted_assemblies = [
            assemblies[i:i + threads] for i in range(0, len(assemblies), threads)
        ]
        for subset in subsetted_assemblies:
            Parallel(n_jobs=threads, prefer="threads")(delayed(run_varifier)(sample,
                                                                                    regions_covered[os.path.basename(sample).replace(".fasta", "")],
                                                                                    input[1],
                                                                                    os.path.join(output[0], os.path.basename(sample).replace(".fasta", "")),
                                                                                    params.container_dir) for sample in subset)

rule assess_assemblies:
    input:
        viridian_assemblies=rules.run_viridian.output,
        simulated_genomes=rules.mask_assemblies.output,
        reference_genome=config["reference_genome"],
        truth_vcf=rules.truth_vcfs.output,
        split_amplicons=rules.split_amplicons.output
    params:
        container_dir=config["container_directory"]
    output:
        directory('varifier_statistics')
    run:
        directories = [output[0], os.path.join(output[0], "ART_assemblies"), os.path.join(output[0], "Badread_assemblies")]
        if not os.path.exists(output[0]):
            for folder in directories:
                os.mkdir(folder)
        art_assemblies = glob.glob(input[0] + '/ART_assemblies/*')
        badread_assemblies = glob.glob(input[0] + '/Badread_assemblies/*')
        # import amplicon statistics
        with open(os.path.join(input[4], 'amplicon_statistics.pickle'), 'rb') as statIn:
            amplicon_stats = pickle.load(statIn)
        # store amplicon regions for each error mode
        dimers = []
        reversions = []
        dropouts = []
        reduced_coverage = []
        # store snp positions and whether it was correctly identifed
        snp_positions = []
        snp_called = []
        # iterate through assemblies and run varifier
        for assem in art_assemblies + badread_assemblies:
            vcf_file = os.path.join(assem, "variants.vcf")
            truth_vcf = os.path.join(input[3], os.path.basename(assem), "04.truth.vcf")
            truth_genome = os.path.join(input[1], os.path.basename(assem) + ".fasta")
            output_dir = os.path.join(output[0], os.path.join(assem.split("/")[-2], os.path.basename(assem)))
            varifier_command = "singularity run " + params.container_dir + "/varifier/varifier.img vcf_eval --truth_vcf " + truth_vcf + " "
           # varifier_command = 'singularity exec "docker://quay.io/iqballab/varifier" varifier vcf_eval --truth_vcf ' + truth_vcf + ' '
            varifier_command += truth_genome + " " + input[2] + " " + vcf_file + " " + output_dir
            shell(varifier_command)
            # import vcf file
            vcf_in = pysam.VariantFile(truth_vcf)
            records = vcf_in.fetch()
            pos = []
            for record in records:
                pos.append(int(record.pos))

               # print(assem, record.ref, record.alts, record.pos)
            # extract samples error modes
            sample_stats = amplicon_stats[os.path.basename(assem)]
            for amplicon in list(sample_stats.keys())[:-1]:
                if sample_stats[amplicon]["has_error"]:
                    if sample_stats[amplicon]["errors"][0] == "primer_dimer":
                        dimers.append((sample_stats[amplicon]["amplicon_start"], sample_stats[amplicon]["amplicon_end"]))
                    if sample_stats[amplicon]["errors"][0] == "primer_reversion":
                        reversions.append((sample_stats[amplicon]["amplicon_start"], sample_stats[amplicon]["amplicon_end"]))
                    if sample_stats[amplicon]["errors"][0] == "random_dropout":
                        dropouts.append((sample_stats[amplicon]["amplicon_start"], sample_stats[amplicon]["amplicon_end"]))
                    if sample_stats[amplicon]["errors"][0] == "primer_SNP":
                        reduced_coverage.append((sample_stats[amplicon]["amplicon_start"], sample_stats[amplicon]["amplicon_end"]))

rule clean_outputs:
    input:
        VGsim_output=rules.VGsim_tree.output,
        phastSim_output=rules.phastSim_evolution.output,
        split_sequences=rules.split_sequences.output,
        split_amplicons=rules.split_amplicons.output,
        read_output=rules.simulate_reads.output,
        cat_reads=rules.cat_reads.output,
        run_viridian=rules.run_viridian.output,
        simulated_genomes=rules.mask_assemblies.output,
        viridian_eval=rules.assess_assemblies.output,
        truth_vcfs=rules.truth_vcfs.output
    shell:
        'rm -rf {input.VGsim_output} {input.phastSim_output} {input.split_sequences} {input.split_amplicons} \
                {input.read_output} {input.cat_reads} {input.run_viridian} {input.simulated_genomes} {input.viridian_eval} \
                {input.truth_vcfs}'
