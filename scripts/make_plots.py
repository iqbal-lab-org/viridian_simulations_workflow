import glob
from joblib import Parallel, delayed
import json
import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import pysam
import seaborn as sns
from statistics import mean
import subprocess
from tqdm import tqdm

def run_cte(primer_scheme,
            assemblies,
            outdir,
            truth_vcf_dir,
            container_dir):
    """run covid-truth-eval"""
    # define the primer scheme
    if primer_scheme == "V3":
        scheme = "COVID-ARTIC-V3"
    if primer_scheme == "V4":
        scheme = "COVID-ARTIC-V4"
    if primer_scheme == "V4.1":
        scheme = "COVID-MIDNIGHT-1200"
    # iterate through assemblies
    manifest = ["name\ttruth_vcf\teval_fasta\tprimers"]
    for assem in assemblies:
        vcf_file = os.path.join(truth_vcf_dir, os.path.basename(assem), "04.truth.vcf")
        manifest.append(os.path.basename(assem) + "\t" + vcf_file + "\t" + os.path.join(assem, "consensus.fa") + "\t" + scheme)
    # save metadata as tsv file
    with open(os.path.join(outdir, "manifest.tsv"), "w") as manifestOut:
        manifestOut.write("\n".join(manifest))
    # run covid-truth-eval
    shell_command = "singularity run " + container_dir + "/covid-truth-eval/cte.img eval_runs \
        --outdir " + os.path.join(outdir, "OUT") + " " + os.path.join(outdir, "manifest.tsv")
    subprocess.run(shell_command, shell=True, check=True)

def run_varifier(assembly,
                truth_vcf_dir,
                reference_genome,
                simulated_genomes,
                output,
                container_dir):
    """Run varifier on viridian assembly to verify SNP calls"""
    vcf_file = os.path.join(assembly, "variants.vcf")
    truth_vcf = os.path.join(truth_vcf_dir, os.path.basename(assembly), "04.truth.vcf")
    truth_genome = os.path.join(simulated_genomes, os.path.basename(assembly) + ".fasta")
    output_dir = os.path.join(output, os.path.join(assembly.split("/")[-2], os.path.basename(assembly)))
    varifier_command = "singularity run " + container_dir + "/varifier/varifier.img vcf_eval --truth_vcf " + truth_vcf + " "
    # varifier_command = 'singularity exec "docker://quay.io/iqballab/varifier" varifier vcf_eval --truth_vcf ' + truth_vcf + ' '
    varifier_command += truth_genome + " " + reference_genome + " " + vcf_file + " " + output_dir
    subprocess.run(varifier_command, shell=True, check=True)

def generate_plots(assembly_list,
                    amplicon_stats,
                    truth_vcf_dir,
                    reference_genome,
                    simulated_genomes,
                    source,
                    threads,
                    output,
                    container_dir):
    """Generate summary plots for the viridian assemblies"""
    # store snp positions and whether it was correctly identifed
    plot_dict = {}
    percentage_dict = {}
    # parallelise varifier
    subsetted_assemblies = [
        assembly_list[i:i + threads] for i in range(0, len(assembly_list), threads)
    ]
    for subset in subsetted_assemblies:
        Parallel(n_jobs=threads, prefer="threads")(delayed(run_varifier)(assem,
                                                                        truth_vcf_dir,
                                                                        reference_genome,
                                                                        simulated_genomes,
                                                                        output,
                                                                        container_dir) for assem in subset)
    # iterate through assemblies
    for assembly in assembly_list:
        output_dir = os.path.join(output, os.path.join(assembly.split("/")[-2], os.path.basename(assembly)))
        # import truth vcf
        truth_vcf_in = pysam.VariantFile(os.path.join(output_dir, "recall", "recall.vcf"))
        truth_records = truth_vcf_in.fetch()
        truth_positions = []
        for record in truth_records:
            truth_positions.append(int(record.pos))
        # import varified vcf
        varifier_vcf_in = pysam.VariantFile(os.path.join(output_dir, "variants_to_eval.filtered.vcf"))
        varifier_records = varifier_vcf_in.fetch()
        varified_positions = []
        for record in varifier_records:
            varified_positions.append(int(record.pos))
        # extract samples error modes
        sample_stats = amplicon_stats[os.path.basename(assembly)]
        amplicon_region = []
        error = []
        # import viridian summary stats to make amplicon starts and ends relative to the reference sequence
        with open(os.path.join(assembly, "log.json")) as inJson:
            summary = json.loads(inJson.read())["read_sampling"]
        for amplicon in list(sample_stats.keys()):
            # make amplicon start and ends relative to ref
            amplicon_region.append((summary[amplicon.split("_LEFT")[0]]["start"]-1, summary[amplicon.split("_LEFT")[0]]["end"]))
            #amplicon_region.append((sample_stats[amplicon]["amplicon_start"], sample_stats[amplicon]["amplicon_end"]))
            if sample_stats[amplicon]["has_error"]:
                error.append(sample_stats[amplicon]["errors"][0])
            else:
                error.append(False)
        regions_seen = []
        # see which varified SNPs are missing
        for position in truth_positions:
            # make position relative to the length of the amplicon
            for region in range(len(amplicon_region)):
                error_mode = error[region]
                if not error_mode:
                    if error[region + 1] == "primer_dimer" or error[region-1] == "primer_dimer":
                        error_mode = "neighbour_primer_dimer"
                    elif error[region + 1] == "primer_SNP" or error[region-1] == "primer_SNP":
                        error_mode = "neighbour_primer_SNP"
                    elif error[region + 1] == "primer_reversion" or error[region-1] == "primer_reversion":
                        error_mode = "neighbour_primer_reversion"
                    elif error[region + 1] == "random_dropout" or error[region-1] == "random_dropout":
                        error_mode = "neighbour_random_dropout"
                    else:
                        error_mode = "no_error"
                if not (position, error_mode) in regions_seen:
                    if position <= amplicon_region[region][1] and position >= amplicon_region[region][0]:
                        # want to look at position calls across all samples
                        # intiate plot dictionary
                        if not error_mode in plot_dict:
                            plot_dict[error_mode] = {}
                            plot_dict[error_mode]["positions"] = {}
                        if not position in plot_dict[error_mode]["positions"]:
                            plot_dict[error_mode]["positions"][position] = {"calls": []}
                        # want to make plots of calls across the amplicon lengths too
                        percentage = round((position - amplicon_region[region][0])*100/(amplicon_region[region][1]-amplicon_region[region][0]), 2)
                        if not error_mode in percentage_dict:
                            percentage_dict[error_mode] = {}
                            percentage_dict[error_mode]["positions"] = {}
                        if not percentage in percentage_dict[error_mode]["positions"]:
                            percentage_dict[error_mode]["positions"][percentage] = {"calls": []}
                        if position in varified_positions:
                            plot_dict[error_mode]["positions"][position]["calls"].append(1)
                            percentage_dict[error_mode]["positions"][percentage]["calls"].append(1)
                        else:
                            plot_dict[error_mode]["positions"][position]["calls"].append(0)
                            percentage_dict[error_mode]["positions"][percentage]["calls"].append(0)
                        regions_seen.append((position, error_mode))
    # add colours to the plot dict
    for to_plot in [plot_dict, percentage_dict]:
        if "primer_dimer" in to_plot:
            to_plot["primer_dimer"]["colour"] = "g"
        if "primer_SNP" in to_plot:  
            to_plot["primer_SNP"]["colour"] = "red"
        if "primer_reversion" in to_plot:
            to_plot["primer_reversion"]["colour"] = "b"
        if "random_dropout" in to_plot:
            to_plot["random_dropout"]["colour"] = "black"
        if "no_error" in to_plot:
            to_plot["no_error"]["colour"] = "purple"
        if "neighbour_primer_dimer" in to_plot:
            to_plot["neighbour_primer_dimer"]["colour"] = "lime"
        if "neighbour_primer_SNP" in to_plot:  
            to_plot["neighbour_primer_SNP"]["colour"] = "magenta"
        if "neighbour_primer_reversion" in to_plot:
            to_plot["neighbour_primer_reversion"]["colour"] = "cyan"
        if "neighbour_random_dropout" in to_plot:
            to_plot["neighbour_random_dropout"]["colour"] = "grey"
    # generate plots
    plt.style.use('seaborn-whitegrid')
    for error in tqdm(plot_dict):
        positions = []
        calls = []
        frequencies = []
        fig = plt.figure()
        ax = fig.add_subplot()
        for site in plot_dict[error]["positions"]:
            positions.append(str(site))
            calls.append(mean(plot_dict[error]["positions"][site]["calls"]))
            frequencies.append(-len(plot_dict[error]["positions"][site]["calls"])/len(assembly_list))
        ordered_positions = positions
        ordered_positions.sort(key=int)
        ordered_positions = [str(x) for x in positions]
        ordered_calls = [calls[positions.index(i)] for i in ordered_positions]
        ax.bar(ordered_positions, ordered_calls, color=plot_dict[error]["colour"])
        # we also are adding the frequency of SNPs at each position below the bars
        ax.bar(ordered_positions, frequencies, color="plum")
        ax.set_xticklabels(ordered_positions, rotation=90, ha='right')
        ax.tick_params(axis='x', which='major', labelsize=5)
        ax.set_ylim([-1.1, 1.1])
        ticks =  ax.get_yticks()
        ax.set_yticklabels([abs(tick) for tick in ticks])
        plt.savefig(os.path.join(output, source + "_" + error + "_SNP_positions.pdf"))

    for error in tqdm(percentage_dict):
        percentages = []
        calls = []
        frequencies = []
        fig = plt.figure()
        ax = fig.add_subplot()
        for pos in percentage_dict[error]["positions"]:
            percentages.append(int(pos))
            calls.append(mean(percentage_dict[error]["positions"][pos]["calls"]))
            frequencies.append(-len(percentage_dict[error]["positions"][pos]["calls"])/len(assembly_list))
        #ordered_percentages = natsorted(percentages)
        #ordered_calls = [calls[percentages.index(i)] for i in ordered_percentages]
        ax.scatter(percentages, calls, s=10, color=percentage_dict[error]["colour"])
        # we also are adding the frequency of SNPs at each position below the bars
        ax.bar(percentages, frequencies, color="plum")
        ax.set_xlim([-0.1, 100.1])
        ax.set_ylim([-1.1, 1.1])
        ticks =  ax.get_yticks()
        ax.set_yticklabels([abs(tick) for tick in ticks])
        plt.savefig(os.path.join(output, source + "_" + error + "_percentage_SNP_positions.pdf"))

def generate_heatmap(eval_dir, method, output_dir):
    """Generate heatmap of viridian assembly SNPs vs ref"""
    cte_files = glob.glob(os.path.join(eval_dir, "*"))
    print(cte_files)
    call_dict = {}
    for subdir in tqdm(cte_files):
        per_position = pd.read_csv(os.path.join(subdir, "per_position.tsv"), sep='\t')
        for pos in range(len(per_position["Ref_pos"])):
            if not per_position["Ref_pos"][pos] in call_dict:
                call_dict[per_position["Ref_pos"][pos]] = []
            else:
                #if any(per_position["Consensus_category"][pos] == cat for cat in ["Called_ref", "Called_correct_alt", "Called_correct_IUPAC"]):
                call = 1
                if any(per_position["Consensus_category"][pos] == cat for cat in ["Called_wrong_alt", "Called_wrong_IUPAC", "Called_wrong_indel", "Called_N"]):
                    call = 0
                call_dict[per_position["Ref_pos"][pos]].append(call)
    heat_values = []
    positions = []
    for pos in call_dict:
        if not mean(call_dict[pos]) == 1:
            heat_values.append(mean(call_dict[pos]))
            positions.append(pos)
    #plt.imshow(np.array(heat_values).reshape(1, len(heat_values)), cmap='plasma')
    heat_matrix = np.matrix(heat_values)
    fig, ax = plt.subplots()
    #hm = sns.heatmap(heat_matrix, cbar=True, cmap="plasma", fmt='.2f', square=True, xticklabels=False, yticklabels=False)
    plt.plot(positions, heat_values)
    plt.savefig(os.path.join(output_dir, method + "_covid_truth_eval_heatmap.pdf"))
