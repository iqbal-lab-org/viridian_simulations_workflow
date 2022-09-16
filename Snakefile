import glob
import gzip
import os
import pickle
import pysam
import random
import shutil
from tqdm import tqdm

from scripts.error_modes import clean_genome, find_primer_scheme
from scripts.make_plots import run_cte, run_varifier

configfile: 'config.yml'

def list_samples():
    return [os.path.basename(s) for s in glob.glob(os.path.join("amplicon_sequences", "*")) if not "/amplicon" in s and not "tmp" in s]

def viridian_art_input(wildcards):
    checkpoint_output = checkpoints.split_amplicons.get(**wildcards).output[0]
    return expand("viridian_ART_assemblies/{SAMPLE}", SAMPLE=list_samples())

def viridian_badread_input(wildcards):
    checkpoint_output = checkpoints.split_amplicons.get(**wildcards).output[0]
    return expand("viridian_Badread_assemblies/{SAMPLE}", SAMPLE=list_samples())

def artic_art_input(wildcards):
    checkpoint_output = checkpoints.split_amplicons.get(**wildcards).output[0]
    return expand("artic_ART_assemblies/{SAMPLE}", SAMPLE=list_samples())

def artic_badread_input(wildcards):
    checkpoint_output = checkpoints.split_amplicons.get(**wildcards).output[0]
    return expand("epi2me_Badread_assemblies/{SAMPLE}", SAMPLE=list_samples())

def truth_vcf_input(wildcards):
    checkpoint_output = checkpoints.split_amplicons.get(**wildcards).output[0]
    return expand("truth_vcfs/{SAMPLE}", SAMPLE=list_samples())

rule all:
    input:
        viridian_art_input,
        viridian_badread_input,
        artic_art_input,
        artic_badread_input,
        truth_vcf_input,
        "cte_viridian_output",
        #"cte_artic_output",
        "usher_phylogenies/artic_ART_phylogeny",
        #"usher_phylogenies/epi2me_Badread_phylogeny",
        "usher_phylogenies/viridian_ART_phylogeny",
        "usher_phylogenies/viridian_Badread_phylogeny"

rule phastSim_evolution:
    input:
        tree_file="reformatted_tree.nwk",
        reference_genome=config["reference_genome"]
    output:
        directory(config['phastSim']['output_dir'])
    params:
        seed=config['seed'],
        container_dir=config["container_directory"],
        rate_parameter=config['phastSim']["rate_parameter"]
    shell:
        'mkdir {output} && singularity run {params.container_dir}/images/phastSim.img --outpath {output}/ --seed {params.seed} --createFasta \
            --createInfo --createNewick --createPhylip --scale 0.0001 --treeFile {input.tree_file} --hyperMutProbs 0.001 0.001 --hyperMutRates 2.0 5.0 \
            --invariable 0.1 --alpha {params.rate_parameter} --omegaAlpha {params.rate_parameter} --codon \
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

checkpoint split_amplicons:
    input:
        split_sequences=rules.split_sequences.output,
        phastSim_dir=rules.phastSim_evolution.output
    output:
        directory('amplicon_sequences')
    params:
        primer_scheme=config['primer_scheme'],
        scheme_dir=config["scheme_dir"],
        seed=config["seed"],
        random_dropout_probability=config['split_amplicons']['random_dropout_probability'],
        primer_dimer_probability=config['split_amplicons']['primer_dimer_probability'],
        match_coverage_mean=config['split_amplicons']['match_coverage_mean'],
        match_coverage_sd=config['split_amplicons']['match_coverage_sd'],
        mismatch_coverage_mean=config['split_amplicons']['mismatch_coverage_mean'],
        mismatch_coverage_sd=config['split_amplicons']['mismatch_coverage_sd'],
        container_dir=config["container_directory"]
    threads: config["threads"]
    resources:
	    mem_mb=lambda wildcards, attempt: 1000 * attempt, threads=config["threads"]
    shell:
        "venv/bin/python scripts/error_modes.py --scheme {params.primer_scheme} --scheme-dir {params.scheme_dir} --input-dir {input.split_sequences} \
            --output-dir {output} --container-dir {params.container_dir} --seed {params.seed} --dropout-prob {params.random_dropout_probability} \
            --dimer-prob {params.primer_dimer_probability} --match-mean {params.match_coverage_mean} --match-sd {params.match_coverage_sd} \
            --mismatch-mean {params.mismatch_coverage_mean} --mismatch-sd {params.mismatch_coverage_sd} --threads {threads}"

rule mask_assemblies:
    input:
        amplicon_sequences="amplicon_sequences/{SAMPLE}",
        simulated_genomes="simulated_genomes"
    output:
        "masked_truth_assemblies/{SAMPLE}.fasta"
    params:
        mask_assemblies=config['mask_assemblies']["apply_mask"]
    threads: 1
    run:
        # check the input sample matches the output sample
        assert os.path.basename(input[0]) == os.path.basename(output[0].replace(".fasta", ""))
        # make directory
        if not os.path.exists(os.path.dirname(output[0])):
            os.mkdir(os.path.dirname(output[0]))
        # mask bases outside of the amplicon scheme
        with open(os.path.join(os.path.dirname(input[0]), 'amplicon_statistics.pickle'), 'rb') as statIn:
            amplicon_stats = pickle.load(statIn)
        # mask bases outside of scheme per sample
        sample = os.path.basename(input[0])
        sample_stats = amplicon_stats[sample]
        all_amplicons = list(sample_stats.keys())
        first_amplicon_start = amplicon_stats[sample][all_amplicons[0]]["amplicon_start"]
        last_amplicon_end = amplicon_stats[sample][all_amplicons[-1]]["amplicon_end"]
        # import truth genome and cut off bases at the ends
        sample_name, sample_sequence = clean_genome(os.path.join(input[1], os.path.basename(input[0]) + ".fasta"))
        sample_sequence = list(sample_sequence)
        sample_sequence = sample_sequence[first_amplicon_start: last_amplicon_end]
        # if mask sequences off then just copy the unmasked sequences
        if not params.mask_assemblies:
            pass
        else:
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
        with open(output[0], "w") as outGen:
            outGen.write("\n".join([">" + sample_name, "".join(sample_sequence)]))

rule truth_vcfs:
    input:
        "masked_truth_assemblies/{SAMPLE}.fasta",
        config["reference_genome"],
        "amplicon_sequences"
    output:
        sample_truth=directory("truth_vcfs/{SAMPLE}"),
    params:
        container_dir=config["container_directory"],
        primer_scheme=config['primer_scheme']
    threads: 1
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
            #with open(os.path.join(output_dir, "04.truth.vcf"), "r") as truth_vcf_in:
            with open(os.path.join(output_dir, "04.truth.vcf"), "r") as truth_vcf_in:
                truth_vcf = truth_vcf_in.read()
            to_add = []
            variant_count = int(truth_vcf.splitlines()[-1].split("\t")[2])
            for dropped in dropped_amplicons:
                start_pos = primer_df.loc[primer_df['name'] == dropped[0]].reset_index(drop=True)["ref_start"][0]
                end_pos = primer_df.loc[primer_df['name'] == dropped[1]].reset_index(drop=True)["ref_end"][0]
                variant_count += 1
                variant_line = "MN908947.3\t" + str(start_pos) + "\t" + str(variant_count)
                variant_line += "\tG\tN\t.\tDROPPED_AMP\tAMP_START=" + str(int(start_pos)-1) + ";AMP_END=" + str(int(end_pos)-1) + "\tGT\t1/1\n"
                truth_vcf += variant_line
            truth_vcf = truth_vcf.split("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample")
            truth_vcf_variants = truth_vcf[1].splitlines()[1:]
            first_line = [truth_vcf[1].splitlines()[1]]
            truth_vcf_variants.sort(key=lambda x: int(x.split("\t")[1]))
            truth_vcf[1] = "\n".join(truth_vcf_variants)
            truth_vcf = "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tsample\n".join(truth_vcf) + "\n"
            with open(os.path.join(output_dir, "04.truth_dropped.vcf"), "w") as truth_vcf_out:
                truth_vcf_out.write("\n".join(truth_vcf.splitlines()[1:]))

        # check the truth vcf input sample matches the output
        assert os.path.basename(input[0].replace(".fasta", "")) == os.path.basename(output[0])
        # make directory
        if not os.path.exists(os.path.dirname(output[0])):
            os.mkdir(os.path.dirname(output[0]))
        # import the amplicon statistics file to extract what parts of the assembly are covered by amplicons
        with open(os.path.join(input[2], 'amplicon_statistics.pickle'), 'rb') as statIn:
            amplicon_stats = pickle.load(statIn)
        regions_covered = {}
        amplicons_dropped = {}
        sample = os.path.basename(input[0]).replace(".fasta", "")
        amplicons = list(amplicon_stats[sample].keys())
        regions_covered[sample] = {"start": str(amplicon_stats[sample][amplicons[0]]["amplicon_start"]),
                                "end": str(amplicon_stats[sample][amplicons[len(amplicons)-1]]["amplicon_end"])}
        # record amplicons dropped for each sample so dropped amplicons can be marked in the truth vcf
        amplicons_dropped[sample] = []
        for a in amplicons:
            if amplicon_stats[sample][a]["has_error"] and \
                    any(amplicon_stats[sample][a]["errors"][0] == mode for mode in ["primer_SNP", "random_dropout", "primer_dimer"]):
                amplicons_dropped[sample].append((a.split("---")[0], a.split("---")[1]))
        # import the primer scheme df to get amplicon positions relative to ref
        primer_df, pool1_primers, pool2_primers = find_primer_scheme(params.primer_scheme,
                                                                    "primer_schemes")
        # parallelise make_truth_vcf
        run_varifier(input[0],
                    regions_covered[sample],
                    amplicons_dropped[sample],
                    primer_df,
                    input[1],
                    output[0],
                    params.container_dir)

rule simulate_art_reads:
    input:
        "amplicon_sequences/{SAMPLE}"
    output:
        directory('ART_read_output/{SAMPLE}')
    params:
        read_length=config["simulate_reads"]["illumina_read_length"],
        seed=config["seed"],
        container_dir=config["container_directory"]
    threads: 1
    resources:
	    mem_mb=lambda wildcards, attempt: 1000 * attempt
    log:
        "logs/simulate_art/{SAMPLE}.log"
    run:
        def simulate_ART_reads(genome,
                                output_dir,
                                sample_coverages,
                                read_length,
                                container_dir):
            """Function to run ART on amplicon sequences per simulated genomic sequence"""
            sample_name = os.path.basename(genome)
            if not os.path.exists(output_dir):
                os.mkdir(output_dir)
            for amplicon in sample_coverages[sample_name]:
                try:
                    coverage = sample_coverages[sample_name][amplicon]
                except:
                    continue
                amplicon_file = os.path.join(genome, amplicon + '.fasta')
                read_file = os.path.join(output_dir, amplicon)
                shell_command = 'singularity run ' + container_dir +'/images/ART.img --quiet -amp -p -sam -na -i ' + amplicon_file + \
                        ' -l ' + read_length + ' -f ' + str(coverage) + ' -o ' + read_file
                shell(shell_command)

        # check the input sample matches the output sample
        assert os.path.basename(input[0]) == os.path.basename(output[0])
        # make output dirs
        output_dirs = [os.path.dirname(output[0]), output[0]]
        for o in output_dirs:
            if not os.path.exists(o):
                os.mkdir(o)
        # get rid of undeleted temp files
        shell("rm -rf amplicon_sequences/tmp*")
        # import coverages
        with open(os.path.join(os.path.dirname(input[0]), 'amplicon_coverages.pickle'), 'rb') as coverageHandle:
            sample_coverages = pickle.load(coverageHandle)
        # run ART
        simulate_ART_reads(input[0],
                        output[0],
                        sample_coverages,
                        str(params.read_length),
                        params.container_dir)

rule simulate_badread_reads:
    input:
        "amplicon_sequences/{SAMPLE}"
    output:
        directory('Badread_read_output/{SAMPLE}')
    params:
        read_length=config["simulate_reads"]["illumina_read_length"],
        seed=config["seed"],
        container_dir=config["container_directory"]
    threads: 1
    resources:
	    mem_mb=lambda wildcards, attempt: 1000 * attempt
    log:
        "logs/simulate_badread/{SAMPLE}.log"
    run:
        def simulate_badreads(genome,
                            output_dir,
                            sample_coverages,
                            container_dir):
            sample_name = os.path.basename(genome)
            if not os.path.exists(output_dir):
                os.mkdir(output_dir)
            for amplicon in sample_coverages[sample_name]:
                try:
                    coverage = sample_coverages[sample_name][amplicon]
                except:
                    continue
                amplicon_file = os.path.join(genome, amplicon + '.fasta')
                read_file = os.path.join(output_dir, amplicon) + '.fastq.gz'
                # skip badread if the coverage is 0
                if str(coverage) == "0":
                    continue
                shell_command = 'singularity run ' + container_dir + '/images/Badread.img simulate --identity 94,98.5,3 --reference '
                shell_command += amplicon_file + ' --quantity ' + str(coverage) + 'x |         gzip > ' + read_file
                shell(shell_command)
            return

        # check the input sample matches the output sample
        assert os.path.basename(input[0]) == os.path.basename(output[0])
        # make output dirs
        output_dirs = [os.path.dirname(output[0]), output[0]]
        for o in output_dirs:
            if not os.path.exists(o):
                os.mkdir(o)
        # get rid of undeleted temp files
        shell("rm -rf amplicon_sequences/tmp*")
        # import coverages
        with open(os.path.join(os.path.dirname(input[0]), 'amplicon_coverages.pickle'), 'rb') as coverageHandle:
            sample_coverages = pickle.load(coverageHandle)
        # run Badread
        simulate_badreads(input[0],
                            output[0],
                            sample_coverages,
                            params.container_dir)

rule cat_art_reads:
    input:
        "ART_read_output/{SAMPLE}"
    output:
        fw_read="concatenated_ART_reads/{SAMPLE}_1.fastq",
        rv_read="concatenated_ART_reads/{SAMPLE}_2.fastq"
    resources:
        mem_mb=lambda wildcards, attempt: 1000 * attempt
    log:
        "logs/concatenate_art/{SAMPLE}.log"
    run:
        # check the input sample matches the output sample
        assert os.path.basename(input[0]) == os.path.basename(output[0]).replace("_1.fastq", "") and \
            os.path.basename(input[0]) == os.path.basename(output[1]).replace("_2.fastq", "")
        # make output dirs
        if not os.path.exists(os.path.dirname(output[0])):
            os.mkdir(os.path.dirname(output[0]))
        # cat ART reads
        forward_reads = sorted(glob.glob(os.path.join(input[0], '*1.fq')))
        forward_filename = output[0]
        fw_concat_command = 'cat ' + ' '.join(forward_reads) + ' > ' + forward_filename
        reverse_reads = sorted(glob.glob(os.path.join(input[0], '*2.fq')))
        reverse_filename = output[1]
        rv_concat_command = 'cat ' + ' '.join(reverse_reads) + ' > ' + reverse_filename
        # forward reads
        shell(fw_concat_command)
        # reverse reads
        shell(rv_concat_command)

rule cat_badread_reads:
    input:
        "Badread_read_output/{SAMPLE}"
    output:
        "concatenated_Badread_reads/{SAMPLE}.fastq"
    resources:
        mem_mb=lambda wildcards, attempt: 1000 * attempt
    log:
        "logs/concatenate_badread/{SAMPLE}.log"
    run:
        # check the input sample matches the output sample
        assert os.path.basename(input[0]) == os.path.basename(output[0]).replace(".fastq", "")
        # make output dirs
        if not os.path.exists(os.path.dirname(output[0])):
            os.mkdir(os.path.dirname(output[0]))
        # concat reads
        reads = sorted(glob.glob(os.path.join(input[0], '*.fastq.gz')))
        new_filenames = []
        # gunzip read files for concatenation
        for read_file in reads:
            new_name = read_file.replace(".fastq.gz", ".fastq")
            with gzip.open(read_file, 'rb') as f_in:
                with open(new_name, 'wb') as f_out:
                    shutil.copyfileobj(f_in, f_out)
            new_filenames.append(new_name)
        filename = output[0]
        concat_command = 'cat ' + ' '.join(new_filenames) + ' > ' + filename
        shell(concat_command)

rule viridian_art_assemble:
    input:
        fw_read="concatenated_ART_reads/{SAMPLE}_1.fastq",
        rv_read="concatenated_ART_reads/{SAMPLE}_2.fastq"
    output:
        directory("viridian_ART_assemblies/{SAMPLE}")
    resources:
	    mem_mb=lambda wildcards, attempt: 1000 * attempt
    threads: 1
    log:
        "logs/viridian_assemble/{SAMPLE}.log"
    params:
        viridian_container=config["viridian_assemble"]["viridian_container"]
    run:
        def illumina_viridian_workflow(reference_genome,
                                        fw_read,
                                        rv_read,
                                        output,
                                        viridian_container):
            """Function to run viridian on ART read sets"""
            viridian_command = "singularity run " + viridian_container + " run_one_sample \
                    --tech illumina \
                    --ref_fasta " + reference_genome + " \
                    --reads1 " + fw_read + " \
                    --reads2 " + rv_read + " \
                    --outdir " + output + "/"
            shell(viridian_command)

        # check the input sample matches the output sample
        assert os.path.basename(input[0]).replace("_1.fastq", "") == os.path.basename(output[0]) and \
            os.path.basename(input[1]).replace("_2.fastq", "") == os.path.basename(output[0])
        # make the output directory
        if not os.path.exists(os.path.dirname(output[0])):
            os.mkdir(os.path.dirname(output[0]))
        # run viridian_workflow illumina pipeline
        illumina_viridian_workflow("reference_genome.fasta",
                                    input[0],
                                    input[1],
                                    output[0],
                                    params.viridian_container)


rule viridian_badread_assemble:
    input:
        "concatenated_Badread_reads/{SAMPLE}.fastq"
    output:
        directory("viridian_Badread_assemblies/{SAMPLE}")
    resources:
	    mem_mb=lambda wildcards, attempt: 1000 * attempt
    threads: 1
    log:
        "logs/viridian_assemble/{SAMPLE}.log"
    params:
        viridian_container=config["viridian_assemble"]["viridian_container"]
    run:
        def nanopore_viridian_workflow(reference_genome,
                                        sample,
                                        output,
                                        viridian_container):
            """Function to run viridian on nanopore read sets"""
            viridian_command = "singularity run " + viridian_container + " run_one_sample \
                    --tech ont \
                    --ref_fasta " + reference_genome + " \
                    --reads " + sample + " \
                    --outdir " + output + "/"
            shell(viridian_command)
        # check the input sample matches the output sample
        assert os.path.basename(input[0]).replace(".fastq", "") == os.path.basename(output[0])
        # make the output directory
        if not os.path.exists(os.path.dirname(output[0])):
            os.mkdir(os.path.dirname(output[0]))
        # run viridian_workflow illumina pipeline
        nanopore_viridian_workflow("reference_genome.fasta",
                                    input[0],
                                    output[0],
                                    params.viridian_container)

rule artic_art_assemble:
    input:
        fw_read="concatenated_ART_reads/{SAMPLE}_1.fastq",
        rv_read="concatenated_ART_reads/{SAMPLE}_2.fastq",
        amplicon_sequences="amplicon_sequences"
    output:
        directory("artic_ART_assemblies/{SAMPLE}")
    resources:
	    mem_mb=lambda wildcards, attempt: 1000 * attempt, threads=2
    threads: 2
    log:
        "logs/artic_assemble/{SAMPLE}.log"
    params:
        artic_illumina_container=config["artic_assemble"]["illumina_workflow_container"],
        primer_scheme=config["primer_scheme"],
        main_nf=config["artic_assemble"]["main_nf"],
        scheme_dir=config["artic_assemble"]["scheme_url"],
        nextflow_path=config["nextflow_path"]
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
            if not os.path.exists(forward_rd + ".gz"):
                shell("gzip " + forward_rd + " " + reverse_rd)
            forward_rd = forward_rd + ".gz"
            reverse_rd = reverse_rd + ".gz"
            shell_command = "venv/bin/python scripts/run_connor_pipeline.py --sif " + sif_file + " "
            shell_command += "--main_nf " + main_nf + " --outdir " + output_dir + " "
            shell_command += "--scheme_url " + scheme_url + " --scheme_version " + primer_scheme + " "
            shell_command += "--ilm1 " + forward_rd + " --ilm2 " + reverse_rd + " --nextflow_path " + nextflow_path + " "
            shell_command += "--sample_name " + os.path.basename(output_dir)
            if not os.path.exists(forward_rd.replace(".gz", "")):
                shell_command += " && gunzip " + forward_rd + " " + reverse_rd
            shell(shell_command)
            # cut off ends of assemblies that lie outside of the amplicon scheme
            filename = os.path.join(output_dir, "consensus.fa")
            sample_name, sample_sequence = clean_genome(filename)
            sample_sequence = list(sample_sequence)
            sample_sequence = sample_sequence[regions_covered["scheme_start"]: regions_covered["scheme_end"]]
            # write out the masked simulated sequence
            with open(filename.replace("consensus", "consensus_trimmed"), "w") as outGen:
                outGen.write("\n".join([">" + sample_name, "".join(sample_sequence)]))
        # check the input sample matches the output sample
        assert os.path.basename(input[0]).replace("_1.fastq", "") == os.path.basename(output[0]) \
            and os.path.basename(input[1]).replace("_2.fastq", "") == os.path.basename(output[0])
        # make output directory
        if not os.path.exists(os.path.dirname(output[0])):
            os.mkdir(os.path.dirname(output[0]))
        # import the primer scheme df to get amplicon scheme start and end relative to ref
        primer_df, pool1_primers, pool2_primers = find_primer_scheme(params.primer_scheme,
                                                                    "primer_schemes")
        with open(os.path.join(input[2], 'amplicon_statistics.pickle'), 'rb') as statIn:
            amplicon_stats = pickle.load(statIn)
        regions_covered = {}
        sample_stats = amplicon_stats[os.path.basename(input[0]).replace("_1.fastq", "").replace(".gz", "")]
        amplicons = list(sample_stats.keys())
        scheme_start = int(primer_df.loc[primer_df['name'] == amplicons[0].split("---")[0]].reset_index(drop=True)["ref_start"][0])
        scheme_end = int(primer_df.loc[primer_df['name'] == amplicons[len(amplicons)-1].split("---")[1]].reset_index(drop=True)["ref_start"][0]) + 1
        regions_covered = {"scheme_start": scheme_start,
                            "scheme_end": scheme_end}
        del amplicon_stats
        # run artic assembly pipeline
        illumina_artic_assemble(input[0],
                                input[1],
                                params.artic_illumina_container,
                                params.main_nf,
                                params.scheme_dir,
                                output[0],
                                params.nextflow_path,
                                params.primer_scheme,
                                regions_covered)

rule epi2me_badread_assemble:
    input:
        "concatenated_Badread_reads/{SAMPLE}.fastq",
        "amplicon_sequences"
    output:
        directory("epi2me_Badread_assemblies/{SAMPLE}")
    resources:
	    mem_mb=lambda wildcards, attempt: 1000 * attempt
    threads: 1
    log:
        "logs/epi2me/{SAMPLE}.log"
    params:
        main_nf=config["epi2me_badread_assemble"]["main_nf"],
        cache_dir=config["epi2me_badread_assemble"]["nxf_sing_cache"],
        primer_scheme=config["primer_scheme"],
        nextflow_path=config["nextflow_path"]
    run:
        def epi2me_artic_assemble(main_nf,
                                output_dir,
                                scheme,
                                read_file,
                                regions_covered,
                                nextflow_path,
                                nextflow_cache):
            """run nanopore artic nextflow pipeline"""
            if not os.path.exists(read_file + ".gz"):
                shell("gzip " + read_file)
            read_file = read_file + ".gz"
            shell_command = "venv/bin/python scripts/run_epi2me.py --main_nf " + main_nf + " "
            shell_command += "--force --sample_name " + os.path.basename(output_dir) + " "
            shell_command += "--work_root_dir " + output_dir + " --outdir " + output_dir + " "
            shell_command += "--nxf_sing_cache " + nextflow_cache + " "
            shell_command += "--scheme_version ARTIC/" + scheme + " "
            shell_command += "--reads " + read_file + " --nextflow_path " + nextflow_path
            if not os.path.exists(read_file.replace(".gz", "")):
                shell_command += " && gunzip " + read_file
            shell(shell_command)
            # cut off ends of assemblies that lie outside of the amplicon scheme
            filename = os.path.join(output_dir, "consensus.fa")
            sample_name, sample_sequence = clean_genome(filename)
            sample_sequence = list(sample_sequence)
            sample_sequence = sample_sequence[regions_covered["scheme_start"]: regions_covered["scheme_end"]]
            # write out the masked simulated sequence
            with open(filename.replace("consensus", "consensus_trimmed"), "w") as outGen:
                outGen.write("\n".join([">" + sample_name, "".join(sample_sequence)]))
        # check the input sample matches the output sample
        assert os.path.basename(input[0]).replace(".fastq", "") == os.path.basename(output[0])
        # make output directory
        if not os.path.exists(os.path.dirname(output[0])):
            os.mkdir(os.path.dirname(output[0]))
        # import the primer scheme df to get amplicon scheme start and end relative to ref
        primer_df, pool1_primers, pool2_primers = find_primer_scheme(params.primer_scheme,
                                                                    "primer_schemes")
        with open(os.path.join(input[1], 'amplicon_statistics.pickle'), 'rb') as statIn:
            amplicon_stats = pickle.load(statIn)
        regions_covered = {}
        sample_stats = amplicon_stats[os.path.basename(input[0]).replace(".fastq", "").replace(".gz", "")]
        amplicons = list(sample_stats.keys())
        scheme_start = int(primer_df.loc[primer_df['name'] == amplicons[0].split("---")[0]].reset_index(drop=True)["ref_start"][0])
        scheme_end = int(primer_df.loc[primer_df['name'] == amplicons[len(amplicons)-1].split("---")[1]].reset_index(drop=True)["ref_start"][0]) + 1
        regions_covered = {"scheme_start": scheme_start,
                            "scheme_end": scheme_end}
        del amplicon_stats
        # run artic assembly pipeline
        epi2me_artic_assemble(params.main_nf,
                            output[0],
                            params.primer_scheme,
                            input[0],
                            regions_covered,
                            params.nextflow_path,
                            params.cache_dir)


def aggregated_va_assemblies(wildcards):
    checkpoint_output = checkpoints.split_amplicons.get(**wildcards).output[0]
    return expand("viridian_ART_assemblies/{sample}", \
        sample=glob_wildcards(os.path.join("viridian_ART_assemblies", "{sample}")).sample)

def aggregated_vb_assemblies(wildcards):
    checkpoint_output = checkpoints.split_amplicons.get(**wildcards).output[0]
    return expand("viridian_Badread_assemblies/{sample}", \
        sample=glob_wildcards(os.path.join("viridian_Badread_assemblies", "{sample}")).sample)

def aggregated_aa_assemblies(wildcards):
    checkpoint_output = checkpoints.split_amplicons.get(**wildcards).output[0]
    return expand("artic_ART_assemblies/{sample}", \
        sample=glob_wildcards(os.path.join("artic_ART_assemblies", "{sample}")).sample)

def aggregated_eb_assemblies(wildcards):
    checkpoint_output = checkpoints.split_amplicons.get(**wildcards).output[0]
    return expand("epi2me_Badread_assemblies/{sample}", \
        sample=glob_wildcards(os.path.join("epi2me_Badread_assemblies", "{sample}")).sample)

def aggregated_tvs(wildcards):
    checkpoint_output = checkpoints.split_amplicons.get(**wildcards).output[0]
    return expand("truth_vcfs/{sample}", \
        sample=glob_wildcards(os.path.join("truth_vcfs", "{sample}")).sample)

rule viridian_covid_truth_eval:
    input:
        aggregated_va_assemblies,
        aggregated_vb_assemblies,
        aggregated_tvs
    output:
        directory("cte_viridian_output")
    threads: 1
    resources:
	    mem_mb=lambda wildcards, attempt: 1000 * attempt
    params:
        container_dir=config["container_directory"],
        primer_scheme=config["primer_scheme"],
    run:
        # make output dirs
        for sub_dir in [output[0], os.path.join(output[0], "ART_assemblies"), os.path.join(output[0], "Badread_assemblies")]:
            if not os.path.exists(sub_dir):
                os.mkdir(sub_dir)
        # list viridian assemblies
        art_assemblies = sorted([f for f in input if "viridian_ART" in f and len(f.split("/")) == 2])
        badread_assemblies = sorted([f for f in input if "viridian_Badread" in f and len(f.split("/")) == 2])
        truth_vcfs = sorted([f for f in input if "truth_vcfs" in f])
        # run covid truth eval
        run_cte(params.primer_scheme,
                art_assemblies,
                os.path.join(output[0], "ART_assemblies"),
                os.path.dirname(truth_vcfs[0]),
                params.container_dir,
                "viridian")
        run_cte(params.primer_scheme,
                badread_assemblies,
                os.path.join(output[0], "Badread_assemblies"),
                os.path.dirname(truth_vcfs[0]),
                params.container_dir,
                "viridian")

rule artic_covid_truth_eval:
    input:
        aggregated_aa_assemblies,
        aggregated_eb_assemblies,
        aggregated_tvs
    output:
        directory("cte_artic_output")
    threads:
        1
    params:
        container_dir=config["container_directory"],
        primer_scheme=config["primer_scheme"]
    resources:
	    mem_mb=lambda wildcards, attempt: 1000 * attempt
    run:
        # make output dirs
        for sub_dir in [output[0], os.path.join(output[0], "ART_assemblies"), os.path.join(output[0], "Badread_assemblies")]:
            if not os.path.exists(sub_dir):
                os.mkdir(sub_dir)
        # list artic assemblies
        art_assemblies = sorted([f for f in input if "artic_ART" in f and len(f.split("/")) == 2])
        badread_assemblies = sorted([f for f in input if "epi2me_Badread" in f and len(f.split("/")) == 2])
        truth_vcfs = sorted([f for f in input if "truth_vcfs" in f])
        # run covid truth eval
        run_cte(params.primer_scheme,
                art_assemblies,
                os.path.join(output[0], "ART_assemblies"),
                os.path.dirname(truth_vcfs[0]),
                params.container_dir,
                "artic")
        run_cte(params.primer_scheme,
                badread_assemblies,
                os.path.join(output[0], "Badread_assemblies"),
                os.path.dirname(truth_vcfs[0]),
                params.container_dir,
                "artic")

rule artic_art_phylogeny:
    input:
        aggregated_aa_assemblies
    output:
        directory("usher_phylogenies/artic_ART_phylogeny")
    threads:
        config["threads"]
    resources:
	    mem_mb=lambda wildcards, attempt: 1000 * attempt, cpus=config["threads"]
    run:
        # make output dirs
        if not os.path.exists("usher_phylogenies"):
            os.mkdir("usher_phylogenies")
        # list assemblies
        aa_assemblies = [v for v in input if len(v.split("/")) == 2]
        # write tsv of assemblies for ushonium input
        samples_rows = []
        for a in aa_assemblies:
            name = os.path.basename(a)
            filename = os.path.join(a, "consensus_trimmed.fa")
            samples_rows.append(name + "\t" + filename)
        tsv_path = os.path.join("usher_phylogenies", "artic_ART.tsv")
        with open(tsv_path, "w") as o:
            o.write("\n".join(samples_rows))
        # run ushonium to make the phylogeny
        ushonium_command = 'singularity/ushonium/Scripts/make_jsonl.py --title "artic ART" --cpus ' + str(threads) + ' '
        ushonium_command += tsv_path + " " + output[0]
        shell(ushonium_command)

rule epi2me_badread_phylogeny:
    input:
        aggregated_eb_assemblies
    output:
        directory("usher_phylogenies/epi2me_Badread_phylogeny")
    threads:
        config["threads"]
    resources:
	    mem_mb=lambda wildcards, attempt: 1000 * attempt, cpus=config["threads"]
    run:
        # make output dirs
        if not os.path.exists("usher_phylogenies"):
            os.mkdir("usher_phylogenies")
        # list assemblies
        eb_assemblies = [v for v in input if len(v.split("/")) == 2]
        # write tsv of assemblies for ushonium input
        samples_rows = []
        for a in eb_assemblies:
            name = os.path.basename(a)
            filename = os.path.join(a, "consensus_trimmed.fa")
            samples_rows.append(name + "\t" + filename)
        tsv_path = os.path.join("usher_phylogenies", "epi2me_Badread.tsv")
        with open(tsv_path, "w") as o:
            o.write("\n".join(samples_rows))
        # run ushonium to make the phylogeny
        ushonium_command = 'singularity/ushonium/Scripts/make_jsonl.py --title "epi2me Badread" --cpus ' + str(threads) + ' '
        ushonium_command += tsv_path + " " + output[0]
        shell(ushonium_command)

rule viridian_art_phylogeny:
    input:
        aggregated_va_assemblies
    output:
        directory("usher_phylogenies/viridian_ART_phylogeny")
    threads:
        config["threads"]
    resources:
	    mem_mb=lambda wildcards, attempt: 1000 * attempt, cpus=config["threads"]
    run:
        # make output dirs
        if not os.path.exists("usher_phylogenies"):
            os.mkdir("usher_phylogenies")
        # list assemblies
        va_assemblies = [v for v in input if len(v.split("/")) == 2]
        # write tsv of assemblies for ushonium input
        samples_rows = []
        for a in va_assemblies:
            name = os.path.basename(a)
            filename = os.path.join(a, "consensus.fa")
            samples_rows.append(name + "\t" + filename)
        tsv_path = os.path.join("usher_phylogenies", "viridian_ART.tsv")
        with open(tsv_path, "w") as o:
            o.write("\n".join(samples_rows))
        # run ushonium to make the phylogeny
        ushonium_command = 'singularity/ushonium/Scripts/make_jsonl.py --title "viridian ART" --cpus ' + str(threads) + ' '
        ushonium_command += tsv_path + ' ' + output[0]
        shell(ushonium_command)

rule viridian_badread_phylogeny:
    input:
        aggregated_vb_assemblies
    output:
        directory("usher_phylogenies/viridian_Badread_phylogeny")
    threads:
        config["threads"]
    resources:
	    mem_mb=lambda wildcards, attempt: 1000 * attempt, cpus=config["threads"]
    run:
        # make output dirs
        if not os.path.exists("usher_phylogenies"):
            os.mkdir("usher_phylogenies")
        # list assemblies
        vb_assemblies = [v for v in input if len(v.split("/")) == 2]
        # write tsv of assemblies for ushonium input
        samples_rows = []
        for a in vb_assemblies:
            name = os.path.basename(a)
            filename = os.path.join(a, "consensus.fa")
            samples_rows.append(name + "\t" + filename)
        tsv_path = os.path.join("usher_phylogenies", "viridian_Badread.tsv")
        with open(tsv_path, "w") as o:
            o.write("\n".join(samples_rows))
        # run ushonium to make the phylogeny
        ushonium_command = 'singularity/ushonium/Scripts/make_jsonl.py --title "viridian Badread" --cpus ' + str(threads) + ' '
        ushonium_command += tsv_path + ' ' + output[0]
        shell(ushonium_command)