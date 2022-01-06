import json
import numpy as np
import os
import pandas as pd
import pickle
import random
import shutil
import subprocess
import sys
import tempfile
from tqdm import tqdm

def populate_primer_dfs(left_primers, right_primers):
    """Several primers have alternate entries, this function duplicates rows so all primer combinations are considered"""
    left_basenames = np.array(['_'.join(name.split('_')[:2]) for name in list(left_primers['name'])])
    right_basenames = np.array(['_'.join(name.split('_')[:2]) for name in list(right_primers['name'])])
    # generate 2D array of all left and right primer combinations and get indices of primer pairs
    matches = np.where([[left_basenames[i] == right_basenames[j] for j in range(0, len(right_basenames))] \
            for i in range(0, len(left_basenames))])
    # generate extended dataframe considering all possible alternate primer pairs
    populated_left_df = left_primers.iloc[matches[0]].reset_index(drop=True)
    populated_right_df = right_primers.iloc[matches[1]].reset_index(drop=True)
    return populated_left_df, populated_right_df

def get_primer_positions(primer_metadata,
                         primer_df):
    """Get position of each primer and add to df"""
    primer_start = []
    primer_end = []
    primer_strand = []
    primer_names = []
    for line in primer_metadata:
        tab_split = line.split("\t")
        primer_names.append(tab_split[3])
        primer_start.append(tab_split[1])
        primer_end.append(tab_split[2])
        primer_strand.append(tab_split[5])
    starts = []
    ends = []
    strands = []
    for name in primer_df["name"]:
        index = primer_names.index(name)
        starts.append(int(primer_start[index]))
        ends.append(int(primer_end[index]))
        strands.append(primer_strand[index])
    primer_df["start"] = starts
    primer_df["end"] = ends
    primer_df["strand"] = strands
    # separate primers into separate dfs of left and right primers
    mask = primer_df['name'].str.contains('LEFT')
    left_primers = primer_df[mask]
    right_primers = primer_df[~mask]
    primer_sequences = {}
    # this function ensures the left and right primer dataframes are the same length
    left_primers, right_primers = populate_primer_dfs(left_primers, right_primers)
    for index, row in left_primers.iterrows():
        # there are sometimes alternate primers with a shared name so we're naming based on the primer pairs
        primer_name = row['name'] + '---' + right_primers['name'][index]
        primer_sequences[primer_name] = [row['seq'], right_primers['seq'][index]]
    # ensure primer lengths match
    if len(left_primers) == len(right_primers):
        return left_primers, right_primers, primer_sequences
    else:
        raise ValueError('Left and right primers are different lengths')

def clean_genome(fasta_file):
    """Split input fasta into sample name and sequence"""
    with open(fasta_file, 'r') as g:
        fasta = g.read().splitlines()
    # split fasta into sample name and sequence
    sample_name = fasta[0].split('>')[1]
    sample_sequence = ''.join(fasta[1:])
    return sample_name, sample_sequence

def reverse_complement(primer_sequence):
    """Get reverse complement of sequence to detect mismatches"""
    base_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    complement = ''.join([base_dict[base] for base in primer_sequence])[::-1]
    return complement

def align_primers(genomic_sequence,
                  primer_df,
                  output_dir):
    """Align primers to genome using BWA and return start and end positions"""
    aligned_starts = []
    aligned_ends = []
    for index, row in primer_df.iterrows():
        temp_dir = tempfile.mkdtemp(dir=output_dir)
        primer_fasta = os.path.join(temp_dir, row['name'] + '.fasta')
        # write each primer sequence to fasta file
        with open(primer_fasta, 'w') as primeSeq:
            primeSeq.write('>' + row['name'] + '\n' + row['seq'])
        alignment_file = os.path.join(temp_dir, row['name'] + '.sam')
        bam_file = os.path.join(temp_dir, row['name'] + '.bam')
        bed_file = os.path.join(temp_dir, row['name'] + '.bed')
        # aln primer to genomic sequence using bwa, then convert to sam, bam and bed files
        map_command = 'bwa index ' + genomic_sequence + ' && bwa aln ' + genomic_sequence + ' ' + primer_fasta
        map_command += ' > ' + os.path.join(temp_dir, row['name'] + '.sai') + '; bwa samse ' + genomic_sequence
        map_command += ' ' + os.path.join(temp_dir, row['name'] + '.sai') +  ' ' + primer_fasta
        map_command += ' > ' + alignment_file + ' && samtools view -Sb ' + alignment_file + ' > ' + bam_file
        map_command += ' && bamToBed -i ' + bam_file + ' > ' + bed_file
        subprocess.run(map_command, shell=True, check=True)
        with open(bed_file, 'r') as bedIn:
            bed_content = bedIn.read().split('\t')
        # extract start and end position of aligned primer
        aligned_starts.append(int(bed_content[1]))
        aligned_ends.append(int(bed_content[2]))
        shutil.rmtree(temp_dir)
    primer_df['start'] = aligned_starts
    primer_df['end'] = aligned_ends
    return primer_df

def find_primer_positions(genomic_sequence,
                          left_primers,
                          right_primers,
                          output_dir):
    """Find primer start and stop positions by aligning to genome"""
    aligned_left_primers = align_primers(genomic_sequence,
                                         left_primers,
                                         output_dir)
    aligned_right_primers = align_primers(genomic_sequence,
                                          right_primers,
                                          output_dir)
    return aligned_left_primers, aligned_right_primers

def induce_mismatches(amplicon_sequence,
                      overlapping_PBS,
                      sample_sequence,
                      mismatch_probability):
    """Decide whether to induce SNP in the overlap of the adjacent amplicon's primer site"""
    base_choices = ['A', 'C', 'T', 'G']
    # keep track of if a mismatch has been induced
    mismatch = 0
    mismatch_positions = []
    # get start position of overlapping PBS in the amplicon
    overlap_start = amplicon_sequence.find(overlapping_PBS)
    # induce mismatch in overlapping PBS of first amplicon but not the second
    new_overlap = overlapping_PBS
    for base in range(len(overlapping_PBS)):
        if np.random.binomial(n=1, p=(mismatch_probability)) == 1:
            new_overlap = list(new_overlap)
            # pseudo-randomly decide a base to substitute
            temp_choices = base_choices[:]
            temp_choices.remove(new_overlap[base])
            new_overlap[base] = random.choice(temp_choices)
            new_overlap = "".join(new_overlap)
            mismatch += 1
            try:
                mismatch_positions.append(overlap_start + base)
            except:
                raise ValueError("The primer sequence of the adjancent amplicon was not found")
    # we need to replace the PBS of the second primer with the reference conditional on the mismatch
    if not mismatch == 0:
        # replace old overlap with the new overlap in the amplicon and the sample sequence
        amplicon_sequence = amplicon_sequence.replace(overlapping_PBS, new_overlap)
        sample_sequence = sample_sequence.replace(overlapping_PBS, new_overlap)
    return amplicon_sequence, sample_sequence, mismatch_positions

def split_amplicons(genomic_sequence,
                    left_primers,
                    right_primers,
                    mismatch_probability,
                    output_dir,
                    truth_seqs):
    """Split a single genomic sequence into amplicons using a list of \
        left primer start positions and right primer end positions"""
    sample_name, sample_sequence = clean_genome(genomic_sequence)
    # we need to map primers to the simulated genomes so INDELs are properly considered
    aligned_left_primers, aligned_right_primers = find_primer_positions(genomic_sequence,
                                                                        left_primers,
                                                                        right_primers,
                                                                        output_dir)
    amplicon_sequences = {}
    primer_df_len = len(aligned_left_primers['start'])
    # a dictionary to keep track of the frequency and positions of error modes
    amplicon_stats = {}
    # extract target sequences to detect mismatches and save amplicon sequences
    primer_pairs = {}
    # extract 5'->3' amplicon sequence from sample sequence, including primers
    for position in range(primer_df_len-1):
        primer_id = aligned_left_primers['name'][position] + '---' + aligned_right_primers['name'][position]
        amplicon = sample_sequence[aligned_left_primers['start'][position]: aligned_right_primers['end'][position]]
        # get primer information for the primer binding site of the adjacent amplicon
        overlapping_PBS = sample_sequence[aligned_left_primers['start'][position+1]: aligned_left_primers['end'][position+1]]
        overlapping_primer_id = aligned_left_primers['name'][position+1] + '---' + aligned_right_primers['name'][position+1]
        # introduce SNPs in the overlapping region of one amplicon at a fixed probability and replace the overlapping binding site with the primer sequence
        amplicon, sample_sequence, mismatch_positions = induce_mismatches(amplicon,
                                                                          overlapping_PBS,
                                                                          sample_sequence,
                                                                          mismatch_probability)
        # add error mode to stats dictionary if a mismatch has been induced
        if not primer_id in amplicon_stats:
            amplicon_stats[primer_id] = {"overlap_mismatches": mismatch_positions}
        else:
            amplicon_stats[primer_id]["overlap_mismatches"] = mismatch_positions
        # if an error needs to be made with this amplicon, replace the left hand primer binding site with the reference primer sequence
        if not ("errors" in amplicon_stats[primer_id] and "reference_primer" in amplicon_stats[primer_id]["errors"]):
            amplicon_stats[primer_id].update({"has_error": 0, "errors": []})
        else:
            # replaces primer binding site with the reference
            amplicon = amplicon.replace(amplicon[:aligned_left_primers['end'][position]-aligned_left_primers['start'][position]], \
                truth_seqs[primer_id][0])
        # if there needs to be an error in the adjacent amplicon PBS, keep track of it
        if not mismatch_positions == []:
            amplicon_stats[overlapping_primer_id] = {"has_error": 1}
            amplicon_stats[overlapping_primer_id]["errors"] = ["reference_primer"]
        amplicon_sequences[primer_id] = amplicon
        # also extract primer target sequences to detect mismatches
        primer_targets = [amplicon[:int(aligned_left_primers['length'][position])], \
                reverse_complement(amplicon[-int(aligned_right_primers['length'][position]):])]
        primer_pairs[primer_id] = primer_targets
        # save each amplicon as a separate fasta
        if not os.path.exists(os.path.join(output_dir, sample_name)):
            os.mkdir(os.path.join(output_dir, sample_name))
        with open(os.path.join(output_dir, sample_name, primer_id + '.fasta'), 'w') as outFasta:
            outFasta.write('>' + primer_id + '\n' + amplicon)
    # add final amplicon with no induced mismatches
    final_id = aligned_left_primers['name'][primer_df_len-1] + '---' + aligned_right_primers['name'][primer_df_len-1]
    final_amplicon = sample_sequence[aligned_left_primers['start'][primer_df_len-1]: aligned_right_primers['end'][primer_df_len-1]]
    if final_id in amplicon_stats and "errors" in amplicon_stats[final_id] and "reference_primer" in amplicon_stats[final_id]["errors"]:
        final_amplicon = final_amplicon.replace(final_amplicon, truth_seqs[final_id][0])
    else:
        amplicon_stats.update({final_id: {"errors": [], "has_error": 0}})
    amplicon_sequences[final_id] = final_amplicon
    # extract final target
    final_target = [final_amplicon[:int(aligned_left_primers['length'][position])], \
                reverse_complement(final_amplicon[-int(aligned_right_primers['length'][position]):])]
    primer_pairs[final_id] = final_target
    # save final amplicon as a separate fasta
    for folder in [os.path.join(output_dir, sample_name), "updated_sequences"]:
        if not os.path.exists(folder):
            os.mkdir(os.path.join(folder))
    with open(os.path.join(output_dir, sample_name, final_id + '.fasta'), 'w') as outFasta:
        outFasta.write('>' + final_id + '\n' + final_amplicon)
    # save updated sample sequence with mismatches in primer binding sites
    with open(os.path.join("updated_sequences", sample_name + ".fasta"), 'w') as outGenome:
        outGenome.write(">" + sample_name + "\n" + sample_sequence)
    return {sample_name: [primer_pairs, amplicon_stats]}

def determine_coverage(feature_dict,
                       amplicon_error_stats,
                       random_dropout_probability,
                       match_coverage_mean,
                       match_coverage_sd,
                       mismatch_coverage_mean,
                       mismatch_coverage_sd):
    """Sample amplicon coverages based on presence or absence of primer-binding mismatches \
        and set amplicon coverage to 0 at a fixed probability"""
    amplicon_names = []
    coverages = []
    # determine coverage of each amplicon using a normal distribution
    for primer_pair in feature_dict:
        amplicon_names.append(primer_pair)
        # reduce the coverage if there are mismatches in the primer-binding region
        if feature_dict[primer_pair] == 1:
            coverage = -1
            while coverage < 0:
                coverage = int(round(np.random.normal(loc=float(mismatch_coverage_mean),
                                                scale=float(mismatch_coverage_sd))))
            amplicon_error_stats[primer_pair]['has_error'] = 1
            amplicon_error_stats[primer_pair]["errors"].append("primer_SNP")
            coverages.append(coverage)
        elif feature_dict[primer_pair] == 0:
            coverage = -1
            while coverage < 0:
                coverage = int(round(np.random.normal(loc=float(match_coverage_mean),
                                                scale=float(match_coverage_sd))))
            coverages.append(coverage)
        else:
            raise ValueError("Feature dictionary is formatted incorrectly, \
                    primer pair values must denote mistmatch presence as int(0) or int(1)")
    # we're also simulating random dropout of amplicons
    dropout_array = np.random.binomial(n=1, p=(1-float(random_dropout_probability)), size=[len(feature_dict)])
    dropped_coverages = [a*b for a,b in zip(coverages,dropout_array)]
    # generate dictionary of amplion_name: coverages
    amplicon_coverages = {}
    for name in range(len(amplicon_names)):
        amplicon_coverages[amplicon_names[name]] = dropped_coverages[name]
        amplicon_error_stats[amplicon_names[name]]["coverage"] = int(dropped_coverages[name])
        if dropped_coverages[name] == 0:
            amplicon_error_stats[amplicon_names[name]]['has_error'] = 1
            amplicon_error_stats[amplicon_names[name]]["errors"].append("random_dropout")
    return amplicon_coverages, amplicon_error_stats

def main():
    # path of V3 artic primer metadata
    primer_sequences = "V3/nCoV-2019.tsv"
    primer_positions = "V3/nCoV-2019.primer.bed"
    # import tsv file as pandas dataframe
    primer_df = pd.read_csv(primer_sequences, sep='\t')
    # import bed file as list
    with open(primer_positions, "r") as pos:
        primer_metadata = pos.read().splitlines()
    # combine tsv and bed file data into dfs of left and right primers and returns list of primer sequence pairs
    left_primers, right_primers, truth_seqs = get_primer_positions(primer_metadata, primer_df)
    # split input genomic sequence into amplicons by left primer start and right primer end
    sequence_files = ['reference_genome.fasta']
    output_dir = 'amplicon_sequences'
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)
    sys.stderr.write('\nAligning primers and splitting simulated sequences into amplicons\n')
    result = {}
    for genomic_sequence in tqdm(sequence_files):
        # returns a nested dictionary of all primer targets in the sample
        sample_primer_targets = split_amplicons(genomic_sequence,
                                                left_primers,
                                                right_primers,
                                                0.05,
                                                output_dir,
                                                truth_seqs)
        result.update(sample_primer_targets)
    # separate primer sequences and amplicon stats
    observed_primer_targets = {}
    error_stats = {}
    for sample in result.keys():
        observed_primer_targets[sample] = result[sample][0]
        error_stats[sample] = result[sample][1]
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
        # isolate error stats for the sample
        amplicon_error_stats = error_stats[sample]
        # generate dictionary of amplicon coverages, considering mismatches and random amplicon dropout
        amplicon_coverages, amplicon_error_stats = determine_coverage(primer_snps,
                                                                      amplicon_error_stats,
                                                                      0.01,
                                                                      100,
                                                                      20,
                                                                      5,
                                                                      1)
        sample_coverages[sample] = amplicon_coverages
        error_stats[sample] = amplicon_error_stats
    # Pickle the output
    sys.stderr.write("\nPickling the output\n")
    with open(os.path.join(output_dir, 'amplicon_coverages.pickle'), 'wb') as ampOut:
        # Pickle using the highest protocol available.
        pickle.dump(sample_coverages, ampOut, pickle.HIGHEST_PROTOCOL)
    # save amplicon error statistics
    with open(os.path.join(output_dir, 'amplicon_statistics.json'), 'w') as statOut:
        statOut.write(json.dumps(error_stats))

if __name__== '__main__':
    main()
