import glob
import numpy as np
import os
import pandas as pd
import pickle
import sys
from tqdm import tqdm

def populate_primer_dfs(left_primers, right_primers):
    """Several primers have alternate entries, this function duplicates rows so all primer combinations are considered"""
    left_basenames = np.array(['_'.join(name.split('_')[:2]) for name in list(left_primers['name'])])
    right_basenames = np.array(['_'.join(name.split('_')[:2]) for name in list(right_primers['name'])])
    matches = np.where([[left_basenames[i] == right_basenames[j] for j in range(0, len(right_basenames))] \
            for i in range(0, len(left_basenames))])
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
    with open(fasta_file, 'r') as g:
        fasta = g.read().splitlines()
    # split fasta into sample name and sequence
    sample_name = fasta[0].split('>')[1]
    sample_sequence = ''.join(fasta[1:])
    return sample_name, sample_sequence

def reverse_complement(primer_sequence):
    base_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    complement = ''.join([base_dict[base] for base in primer_sequence])[::-1]
    return complement

def split_amplicons(genome, 
                    left_primers,
                    right_primers,
                    output_dir):
    sample_name, sample_sequence = clean_genome(genome)
    # extract 5'->3' amplicon sequence from sample sequence, including primers
    primer_pairs = {}
    for position in range(len(left_primers['start'])):
        primer_id = left_primers['name'][position] + '---' + right_primers['name'][position]
        amplicon = sample_sequence[left_primers['start'][position]: right_primers['end'][position]]
        # also extract primer target sequences to detect mismatches
        primer_targets = [amplicon[:int(left_primers['length'][position])], \
                reverse_complement(amplicon[-int(right_primers['length'][position]):])]
        primer_pairs[primer_id] = primer_targets
        # save each amplicon as a separate fasta
        if not os.path.exists(os.path.join(output_dir, sample_name)):
            os.mkdir(os.path.join(output_dir, sample_name))
        with open(os.path.join(output_dir, sample_name, primer_id + '.fasta'), 'w') as outFasta:
            outFasta.write('>' + primer_id + '\n' + amplicon)
    return {sample_name: primer_pairs}

def determine_coverage(feature_dict,
                       random_dropout_probability,
                       match_coverage_mean,
                       match_coverage_sd,
                       mismatch_coverage_mean,
                       mismatch_coverage_sd):
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
    return amplicon_coverages

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
    sys.stderr.write('\nSplitting simulated sequences into amplicons\n')
    observed_primer_targets = {}
    for genomic_sequence in tqdm(sequence_files):
        # returns a nested dictionary of all primer targets in the sample
        sample_primer_targets = split_amplicons(genomic_sequence,
                                                left_primers,
                                                right_primers,
                                                output_dir)
        observed_primer_targets.update(sample_primer_targets)
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
                                                0.01,
                                                100,
                                                20,
                                                5,
                                                1)
        sample_coverages[sample] = amplicon_coverages
    # Pickle the output 
    sys.stderr.write("\nPickling the output\n")
    with open(os.path.join(output_dir, 'amplicon_coverages.pickle'), 'wb') as ampOut:
        # Pickle using the highest protocol available.
        pickle.dump(sample_coverages, ampOut, pickle.HIGHEST_PROTOCOL)

if __name__== '__main__':
    main()
