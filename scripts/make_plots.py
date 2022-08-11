import glob
from joblib import Parallel, delayed
import json
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import os
import pandas as pd
import pickle
import pysam
import seaborn as sns
from statistics import mean
import subprocess
import tempfile
from tqdm import tqdm

def run_cte(primer_scheme,
            assemblies,
            outdir,
            truth_vcf_dir,
            container_dir,
            method):
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
    count = 0
    for assem in assemblies:
        vcf_file = os.path.join(truth_vcf_dir, os.path.basename(assem), "04.truth_dropped.vcf")
        #vcf_file = os.path.join(truth_vcf_dir, os.path.basename(assem), "04.truth.vcf")
        if method == "artic":
            assembly_file = os.path.join(assem, "consensus_trimmed.fa")
        else:
            #assembly_file = os.path.join(assem, "consensus.fa")
            assembly_file = os.path.join(assem, "masked.fasta")
        #for a in assemblies:
        #   manifest.append(os.path.basename(a) + "_" + str(count) + "\t" + vcf_file + "\t" + os.path.join(a, "consensus.fa") + "\t" + scheme)
        #  count += 1
        manifest.append(os.path.basename(assem) + "\t" + vcf_file + "\t" + assembly_file + "\t" + scheme)
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
                container_dir,
                artic_vcf_dir=None):
    """Run varifier on viridian assembly to verify SNP calls"""
    if not artic_vcf_dir:
        vcf_file = os.path.join(assembly, "variants.vcf")
        output_dir = os.path.join(output, "viridian_assemblies", os.path.join(assembly.split("/")[-2], os.path.basename(assembly)))
    else:
        vcf_file = os.path.join(artic_vcf_dir, "/".join(assembly.split("/")[-2:]), "04.truth.vcf")
        output_dir = os.path.join(output, "artic_assemblies", os.path.join(assembly.split("/")[-2], os.path.basename(assembly)))
    truth_vcf = os.path.join(truth_vcf_dir, os.path.basename(assembly), "04.truth.vcf")
    truth_genome = os.path.join(simulated_genomes, os.path.basename(assembly) + ".fasta")
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
                    container_dir,
                    artic_vcf_dir=None,
                    viridian_assembly_dir=None):
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
                                                                        container_dir,
                                                                        artic_vcf_dir) for assem in subset)
    if not artic_vcf_dir:
        method_dir = "viridian_assemblies"
    else:
        method_dir = "artic_assemblies"
    # collect all snp calls
    all_snp_calls = {}
    # iterate through assemblies
    for assembly in assembly_list:
        output_dir = os.path.join(output, method_dir, os.path.join(assembly.split("/")[-2], os.path.basename(assembly)))
        # import truth vcf
        truth_vcf_in = pysam.VariantFile(os.path.join(truth_vcf_dir, os.path.basename(assembly), "04.truth.vcf"))
        truth_records = truth_vcf_in.fetch()
        truth_positions = []
        truth_snps = []
        for record in truth_records:
            truth_positions.append(int(record.pos))
            truth_snps.append(record.alts)
        # import varified vcf
        varifier_vcf_in = pysam.VariantFile(os.path.join(output_dir, "variants_to_eval.filtered.vcf"))
        varifier_records = varifier_vcf_in.fetch()
        varified_positions = []
        varified_snps = []
        for record in varifier_records:
            varified_positions.append(int(record.pos))
            varified_snps.append(record.alts)
        # extract samples error modes
        sample_stats = amplicon_stats[os.path.basename(assembly)]
        amplicon_region = []
        error = []
        # import viridian summary stats to make amplicon starts and ends relative to the reference sequence
        if not viridian_assembly_dir:
            with open(os.path.join(assembly, "log.json")) as inJson:
                summary = json.loads(inJson.read())["read_sampling"]
        else:
            with open(os.path.join(viridian_assembly_dir, os.path.basename(assembly), "log.json")) as inJson:
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
        for position in range(len(truth_positions)):
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
                if not (truth_positions[position], error_mode) in regions_seen:
                    if truth_positions[position] <= amplicon_region[region][1] and truth_positions[position] >= amplicon_region[region][0]:
                        # want to look at position calls across all samples
                        # intiate plot dictionary
                        if not error_mode in plot_dict:
                            plot_dict[error_mode] = {}
                            plot_dict[error_mode]["positions"] = {}
                        if not truth_positions[position] in plot_dict[error_mode]["positions"]:
                            plot_dict[error_mode]["positions"][truth_positions[position]] = {"calls": []}
                        # want to make plots of calls across the amplicon lengths too
                        percentage = round((truth_positions[position] - amplicon_region[region][0])*100/(amplicon_region[region][1]-amplicon_region[region][0]), 2)
                        if not error_mode in percentage_dict:
                            percentage_dict[error_mode] = {}
                            percentage_dict[error_mode]["positions"] = {}
                        if not percentage in percentage_dict[error_mode]["positions"]:
                            percentage_dict[error_mode]["positions"][percentage] = {"calls": []}
                        if truth_positions[position] in varified_positions and truth_snps[position][0] == varified_snps[varified_positions.index(truth_positions[position])][0]:
                            plot_dict[error_mode]["positions"][truth_positions[position]]["calls"].append(1)
                            percentage_dict[error_mode]["positions"][percentage]["calls"].append(1)
                        else:
                            plot_dict[error_mode]["positions"][truth_positions[position]]["calls"].append(0)
                            percentage_dict[error_mode]["positions"][percentage]["calls"].append(0)
                        regions_seen.append((truth_positions[position], error_mode))
        all_snp_calls[os.path.basename(assembly)] = {"truth_snp_positions": truth_positions, "truth_snp_calls": truth_snps, "varifier_snp_positions": varified_positions, "varifier_snp_calls": varified_snps}

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
        position_call_dicts = []
        position_frequency_dicts = []
        fig = plt.figure()
        ax = fig.add_subplot()
        for site in plot_dict[error]["positions"]:
            position_call_dicts.append({str(site): mean(plot_dict[error]["positions"][site]["calls"])})
            position_frequency_dicts.append({str(site): -len(plot_dict[error]["positions"][site]["calls"])/len(assembly_list)})
        position_call_dicts = sorted(position_call_dicts, key=lambda x: int(list(x.keys())[0])) 
        position_frequency_dicts = sorted(position_frequency_dicts, key=lambda x: int(list(x.keys())[0])) 
        ax.bar([list(position.keys())[0] for position in position_call_dicts], [list(call.values())[0] for call in position_call_dicts], color=plot_dict[error]["colour"])
        # we also are adding the frequency of SNPs at each position below the bars
        ax.bar([list(position.keys())[0] for position in position_call_dicts], [list(freq.values())[0] for freq in position_frequency_dicts], color="plum")
        ax.set_xticklabels([list(position.keys())[0] for position in position_call_dicts], rotation=90, ha='right')
        ax.tick_params(axis='x', which='major', labelsize=5)
        ax.set_ylim([-1.1, 1.1])
        ticks =  ax.get_yticks()
        ax.set_yticklabels([abs(tick) for tick in ticks])
        plt.savefig(os.path.join(output, source + "_" + error + "_SNP_positions.png"), dpi=300)

    for error in tqdm(percentage_dict):
        percentage_call_dicts = []
        percentage_frequency_dicts = []
        fig = plt.figure()
        ax = fig.add_subplot()
        for pos in percentage_dict[error]["positions"]:
            percentage_call_dicts.append({int(pos): mean(percentage_dict[error]["positions"][pos]["calls"])})
            percentage_frequency_dicts.append({int(pos): -len(percentage_dict[error]["positions"][pos]["calls"])/len(assembly_list)})
        ax.scatter([list(percentage.keys())[0] for percentage in percentage_call_dicts], [list(call.values())[0] for call in percentage_call_dicts], s=10, color=percentage_dict[error]["colour"])
        # we also are adding the frequency of SNPs at each position below the bars
        ax.bar([list(percentage.keys())[0] for percentage in percentage_call_dicts], [list(freq.values())[0] for freq in percentage_frequency_dicts], color="plum")
        ax.set_xlim([-0.1, 100.1])
        ax.set_ylim([-1.1, 1.1])
        ticks =  ax.get_yticks()
        ax.set_yticklabels([abs(tick) for tick in ticks])
        plt.savefig(os.path.join(output, source + "_" + error + "_percentage_SNP_positions.png"), dpi=300)
    return all_snp_calls

def generate_heatmap(eval_dir, method, output_dir):
    """Generate heatmap of viridian assembly SNPs vs ref"""
    cte_files = glob.glob(os.path.join(eval_dir, "*"))
    call_dict = {}
    for subdir in tqdm(cte_files):
        per_position = pd.read_csv(os.path.join(subdir, "per_position.tsv"), sep='\t')
      # for category in range(len(per_position["Truth_category"])):
           # if per_position["SNP_true_alt"][category] == "SNP_true_alt"
        for pos in range(len(per_position["Ref_pos"])):
            if not per_position["Ref_pos"][pos] in call_dict:
                call_dict[per_position["Ref_pos"][pos]] = []
            else:
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
    plt.savefig(os.path.join(output_dir, method + "_covid_truth_eval_heatmap.png"), dpi=300)

def plot_varifier_calls(viridian_art_snps,
                        viridian_badread_snps,
                        artic_art_snps,
                        artic_badread_snps,
                        output_dir,):
    """Plot the union of snp calls"""
    # iterate through call set
    all_sets = [viridian_art_snps, viridian_badread_snps, artic_art_snps, artic_badread_snps]
    all_call_dict = {}
    all_calls_union = []
    for calls in range(len(all_sets)):
        wrong = set()
        call_source = all_sets[calls]
        call_by_base = {}
        # iterate through samples
        for sample in call_source:
            calls_union = list(set(call_source[sample]["truth_snp_positions"]) | set(call_source[sample]["varifier_snp_positions"]))
            calls_union.sort(key=int)
            all_calls_union += calls_union
            if not sample in call_by_base:
                call_by_base[sample] = {}
            for position in calls_union:
                if not position in call_by_base[sample]:
                    call_by_base[sample][position] = ""
                if position in call_source[sample]["truth_snp_positions"] and position in call_source[sample]["varifier_snp_positions"]:
                    if call_source[sample]["varifier_snp_calls"][call_source[sample]["varifier_snp_positions"].index(position)][0] == call_source[sample]["truth_snp_calls"][call_source[sample]["truth_snp_positions"].index(position)][0]:
                        call_by_base[sample][position] = 1
                    else:
                        call_by_base[sample][position] = -1
                        wrong.add(position)
                else:
                    call_by_base[sample][position] = -1
                    wrong.add(position)
        if calls == 0:
            all_call_dict["viridian_ART_heatmap.png"] = call_by_base
        elif calls == 1:
            all_call_dict["viridian_Badread_heatmap.png"] = call_by_base
        elif calls == 2:
            all_call_dict["artic_ART_heatmap.png"] = call_by_base
        else:
            all_call_dict["artic_Badread_heatmap.png"] = call_by_base
    # make sure all samples have a value for all bases
    all_calls_union = list(set(all_calls_union))
    all_calls_union.sort(key=int)
    for key in all_call_dict:
        heat_matrix = []
        for samp in all_call_dict[key]:
            heat_values = []
            for pos in all_calls_union:
                if not pos in all_call_dict[key][samp].keys():
                    heat_values.append(0)
                else:
                    heat_values.append(all_call_dict[key][samp][pos])
            heat_matrix.append(heat_values)
        # plot the values
        heat_matrix = np.matrix(heat_matrix)
        fig, ax = plt.subplots()
        myColors = ("r", "w", "b")
        cmap = LinearSegmentedColormap.from_list('Custom', myColors, len(myColors))
        hm = sns.heatmap(heat_matrix, cbar=False, cmap=cmap, vmin=-0.5, vmax=0.5, fmt='.2f', square=True, xticklabels=False, yticklabels=False)
        plt.tight_layout()
        plt.savefig(os.path.join(output_dir, key), dpi=300)

def pairwise_compare(first_assemblies,
                     second_assemblies_dir,
                     output_dir,
                     comparison,
                     container_dir):
    """Get pairwise differences in bases"""
    temp_dir = tempfile.mkdtemp(dir=output_dir)
    alignment_files = []
    for assembly in first_assemblies:
        if comparison == "viridian_artic":
            first_assembly = os.path.join(assembly, "consensus.fa")
            second_assembly = os.path.join(second_assemblies_dir, os.path.basename(assembly), "consensus_trimmed.fa")
        if comparison == "viridian_simulated":
            first_assembly = os.path.join(assembly, "consensus.fa")
            second_assembly = os.path.join(second_assemblies_dir, os.path.basename(assembly) + ".fasta")
        if comparison == "artic_simulated":
            first_assembly = os.path.join(assembly, "consensus_trimmed.fa")
            second_assembly = os.path.join(second_assemblies_dir, os.path.basename(assembly) + ".fasta")
        input_alignment = os.path.join(temp_dir, "unaligned_" + os.path.basename(assembly) + ".fasta")
        sequences = []
        with open(first_assembly) as inFirst:
            first_assembly = inFirst.read().splitlines()[1:]
        with open(second_assembly) as inSecond:
            second_assembly = inSecond.read().splitlines()[1:]
        sequences.append(">" + os.path.basename(assembly) + "_" + comparison.split("_")[0] + "\n" + "".join(first_assembly))
        sequences.append(">" + os.path.basename(assembly) + "_" + comparison.split("_")[1] + "\n" + "".join(second_assembly))
        with open(input_alignment, "w") as outAln:
            outAln.write("\n".join(sequences))
        output_alignment = os.path.join(temp_dir, "aligned_" + os.path.basename(assembly) + ".fasta")
        align_command = container_dir + "/images/Mafft.img --6merpair " + input_alignment + " > " + output_alignment
        subprocess.run(align_command, shell=True, check=True)
        alignment_files.append(output_alignment)
    # compare base by base to build dictionary comparing assembly methods
    base_by_base = {}
    for alignment in alignment_files:
        with open(alignment, "r") as inA:
            content = inA.read()
        methods = content.split(">")[1:]
        for sample in range(len(methods)):
            methods[sample] = "".join(methods[sample].upper().splitlines()[1:])
        first_sequence = list(methods[0])
        second_sequence = list(methods[1])
        for base in range(len(first_sequence)):
            if not str(first_sequence[base] + second_sequence[base]) in base_by_base:
                base_by_base[str(first_sequence[base] + second_sequence[base])] = []
            base_by_base[str(first_sequence[base] + second_sequence[base])].append(base + 1)
    base_counts = {}
    for change in base_by_base:
        base_counts[change] = len(base_by_base[change])
    with open(os.path.join(output_dir, comparison + "_base_by_base.json"), "w") as outJson:
        outJson.write(json.dumps(base_by_base))
    with open(os.path.join(output_dir, comparison + "_base_by_base_counts.json"), "w") as outJson:
        outJson.write(json.dumps(base_counts))

def count_dropped_bases(amplicon_dir,
                        output_dir):
    """Count number of low coverage bases and write out as json"""
    with open(os.path.join(amplicon_dir, 'amplicon_statistics.pickle'), 'rb') as statsHandle:
        amplicon_stats = pickle.load(statsHandle)
    # iterate through samples
    N_stats = {"low_coverage_bases": 0, "number_amplicons": 0}
    for sample in amplicon_stats:
        pos_coverage_dict = {}
        # iterate through amplicons
        for amplicon in amplicon_stats[sample]:
            if amplicon_stats[sample][amplicon]["coverage"] < 10:
                N_stats["number_amplicons"] += 1
                current_start = amplicon_stats[sample][amplicon]["amplicon_start"]
                current_end = amplicon_stats[sample][amplicon]["amplicon_end"]
                amplicon_no = int(amplicon.split("_")[1])
                if not amplicon_no == 0:
                    previous_amplicon = "_".join(["SARS-CoV-2", str(amplicon_no - 1), "LEFT---SARS-CoV-2", str(amplicon_no - 1), "RIGHT"])
                    previous_end = amplicon_stats[sample][previous_amplicon]["amplicon_end"]
                else:
                    previous_end = current_start
                if not amplicon_no == 99:
                    next_amplicon = "_".join(["SARS-CoV-2", str(amplicon_no + 1), "LEFT---SARS-CoV-2", str(amplicon_no + 1), "RIGHT"])
                    next_start = amplicon_stats[sample][next_amplicon]["amplicon_start"]
                else:
                    next_start = current_end
                low_coverage = next_start - previous_end   
                N_stats["low_coverage_bases"] = N_stats["low_coverage_bases"] + low_coverage
    # write low coverage base count as json
    with open(os.path.join(output_dir, "low_coverage_base_count.json"), "w") as outJson:
        outJson.write(json.dumps(N_stats))
