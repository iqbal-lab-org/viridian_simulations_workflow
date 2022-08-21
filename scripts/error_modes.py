import argparse
import glob
from joblib import Parallel, delayed
import json
import numpy as np
from numpy.random import RandomState
import os
import pandas as pd
import pickle
import pysam
import re
import shutil
import subprocess
import sys
import tempfile
from tqdm import tqdm

def find_primer_scheme(scheme,
                    scheme_dir):
    """Determine the primer scheme and load the metadata"""
    if scheme == "V3":
        # path of V3 artic primer metadata
        primer_sequences = os.path.join(scheme_dir, scheme, "nCoV-2019.tsv")
        bed_file = os.path.join(scheme_dir, scheme, "nCoV-2019.scheme.bed")
        # import tsv file as pandas dataframe
        primer_df = pd.read_csv(primer_sequences, sep='\t')
        with open(bed_file, "r") as inBed:
            bed_content = inBed.read().splitlines()
        # add primer start and ends to the primer_df
        ref_starts = []
        ref_ends = []
        for name in primer_df["name"]:
            bed_row = next((s for s in bed_content if name in s), None)
            bed_row = bed_row.split("\t")
            ref_starts.append(bed_row[1])
            ref_ends.append(bed_row[2])
        primer_df["ref_start"] = ref_starts
        primer_df["ref_end"] = ref_ends
        # split primer sequences by their pool
        pool1_primers = []
        pool2_primers = []
        for row in range(len(primer_df["pool"])):
            if "_1" in primer_df["pool"][row]:
                pool1_primers.append(primer_df["seq"][row])
            if "_2" in primer_df["pool"][row]:
                pool2_primers.append(primer_df["seq"][row])
        return primer_df, pool1_primers, pool2_primers
    if "V4" in scheme:
        # path of V4 artic primer metadata
        primer_sequences = os.path.join(scheme_dir, scheme, "SARS-CoV-2.primer.bed")
        # import bed file as string
        with open(primer_sequences, "r") as inV4:
            primers = inV4.read().splitlines()
        # keep track of which primers are in pool1 or pool2
        pool1_primers = []
        pool2_primers = []
        # extract primer names and sequences and convert to df
        primer_dict = {"name": [], "seq": [], "pool": [], "ref_start": [], "ref_end": []}
        for row in primers:
            split = row.split("\t")
            primer_dict["name"].append(split[3])
            primer_dict["seq"].append(split[6])
            primer_dict["ref_start"].append(split[1])
            primer_dict["ref_end"].append(split[2])
            if split[4] == "1":
                pool1_primers.append(split[6])
                primer_dict["pool"].append("nCoV-2019_1")
            if split[4] == "2":
                pool2_primers.append(split[6])
                primer_dict["pool"].append("nCoV-2019_2")
        if scheme == "V4":
            return pd.DataFrame(primer_dict), pool1_primers, pool2_primers
        elif scheme == "V4.1":
            with open(os.path.join(scheme_dir, scheme, "SARS-CoV-2.primer_pairs.tsv")) as inTSV:
                pairs = inTSV.read().splitlines()
            primer_data = []
            ordered_primer_names = []
            # the order of the V4.1 primer scheme is a mess
            for pair in pairs:
                split_names = pair.split("\t")
                for name in split_names:
                    name = name.split("_alt")[0]
                    ordered_primer_names += [name, name + "_alt1"]
            # iterate through the correctly ordered names and get the sequence
            primer_pools = []
            for prospective_name in ordered_primer_names:
                if prospective_name in primer_dict["name"]:
                    primer_pools.append(primer_dict["pool"][primer_dict["name"].index(prospective_name)])
                    primer_data.append((prospective_name, primer_dict["seq"][primer_dict["name"].index(prospective_name)], primer_dict["ref_start"][primer_dict["name"].index(prospective_name)], primer_dict["ref_end"][primer_dict["name"].index(prospective_name)]))
            primer_df = pd.DataFrame(primer_data, columns=['name','seq','ref_start','ref_end'])
            primer_df["pool"] = primer_pools
            return primer_df, pool1_primers, pool2_primers
        else:
            ValueError('Only "V3", "V4" and "V4.1" primer schemes are supported')
    else:
        raise ValueError('Only "V3", "V4" and "V4.1" primer schemes are supported')

def correct_alignment(sample_sequence,
                    bed_content,
                    primer_seq):
    """Correct cases where bwa has soft clipped the primer alignment. \
        For now we are assuming there are no INDELs in the primer binding sites"""
    primer_start = int(bed_content[1])
    # get full primer binding sequence
    target_seq = sample_sequence[primer_start: primer_start + len(primer_seq)]
    # compare the primer and binding site to
    mismatches = 0
    mismatch_str = ""
    position_count = 0
    from_end = len(primer_seq)
    for base in range(len(primer_seq)):
        if not primer_seq[base] == target_seq[base]:
            mismatches += 1
            mismatch_str += str(position_count)
            mismatch_str += target_seq[base]
            position_count = 0
            from_end = len(primer_seq) - (base + 1)
        else:
            position_count += 1
    mismatch_str += str(from_end)
    # return a list of tuples
    sam_tags = [("NM", mismatches), ("MD", mismatch_str)]
    # correct end position of primer
    bed_content[2] = primer_start + len(primer_seq)
    return sam_tags, bed_content

def align_primers(genomic_sequence,
                primer_df,
                sample_sequence,
                output_dir,
                container_dir):
    """Align primers to genome using BWA and return start and end positions"""
    aligned_starts = []
    aligned_ends = []
    alignment_stats = {}
    amplicon_dirs = []
    # make temp dir for the sample
    temp_dir = tempfile.mkdtemp(dir=output_dir)
    for index, row in primer_df.iterrows():
        sub_temp_dir = tempfile.mkdtemp(dir=temp_dir)
        primer_fasta = os.path.join(sub_temp_dir, row['name'] + '.fasta')
        # write each primer sequence to fasta file
        with open(primer_fasta, 'w') as primeSeq:
            primeSeq.write('>' + row['name'] + '\n' + row['seq'])
        amplicon_dirs.append(os.path.join(sub_temp_dir, row['name']))
    # write out list of amplicon dirs
    with open(os.path.join(temp_dir, "amplicons.txt"), "w") as AmpDir:
        AmpDir.write("\n".join(amplicon_dirs))
    # parse list of amplicon dirs to reduce singularity overhead
    map_command = "singularity run " + container_dir +"/images/Map.img --dirs " + temp_dir + " --genomic-sequence " + genomic_sequence
    subprocess.run(map_command, shell=True, check=True)
    # iterate through amplicons to get relevant mapped information
    for directory in amplicon_dirs:
        alignment_file = directory + '.sam'
        bed_file = directory + '.bed'
        # open the bed and sam files
        with open(bed_file, 'r') as bedIn:
            bed_content = bedIn.read().split('\t')
        samfile = pysam.AlignmentFile(alignment_file, "r")
        # see if there are mismatches between the primer and the simulated sequence
        for mapped in samfile:
            sam_tags = mapped.tags
            if not mapped.cigar[0][1] == 0:
                sam_tags, bed_content = correct_alignment(sample_sequence,
                                                        bed_content,
                                                        mapped.query_sequence)
            mismatches = [tag[1] for tag in sam_tags if tag[0] == "NM"][0]
            mismatch_dict = {}
            if not mismatches == 0:
                # if there are primer mismatches then get the positions
                positions = [tag[1] for tag in sam_tags if tag[0] == "MD"][0]
                positions = [pos for pos in re.split(r"([A-Z]+)", positions)][:-1]
                absolute_position = 0
                for pos in range(len(positions)-1):
                    if positions[pos].isnumeric():
                        current = int(positions[pos])
                        adjacent_char = positions[pos+1]
                        # need to adjust the .sam positions
                        if not pos == 0:
                            current = absolute_position + current + 1
                        absolute_position += current
                        mismatch_dict[current] = adjacent_char
            alignment_stats[os.path.basename(directory)] = {"mismatches": mismatches,
                                                            "positions": mismatch_dict}
        # extract start and end position of aligned primer
        aligned_starts.append(int(bed_content[1]))
        aligned_ends.append(int(bed_content[2]))
    # store the mapped stop and start positions
    primer_df['start'] = aligned_starts
    primer_df['end'] = aligned_ends
    # remove the temp_dir
    try:
        shutil.rmtree(temp_dir, ignore_errors=True)
    except FileNotFoundError:
        pass
    return primer_df, alignment_stats

def clean_genome(fasta_file):
    """Split input fasta into sample name and sequence"""
    with open(fasta_file, 'r') as g:
        fasta = g.read().splitlines()
    # split fasta into sample name and sequence
    sample_name = os.path.basename(fasta_file).split(".fasta")[0]
    sample_sequence = ''.join(fasta[1:])
    return sample_name, sample_sequence

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

def refine_primers(primer_df,
                alignment_stats):
    """This function selects the best matching left and right primer for each amplicon \
        if both match, then the non-alternative primer is selected"""
    rows_to_drop = []
    for name in range(len(primer_df["name"])-1):
        if "alt" in primer_df["name"][name+1] and (primer_df["start"][name+1] < primer_df["start"][name] + len(primer_df["seq"][name]) \
            or primer_df["end"][name+1] - len(primer_df["seq"][name+1]) > primer_df["end"][name]):
            if alignment_stats[primer_df["name"][name]]["mismatches"] <= alignment_stats[primer_df["name"][name + 1]]["mismatches"]:
                # drop the alternate primer if it is not the best match
                rows_to_drop.append(name+1)
            else:
                rows_to_drop.append(name)
    primer_df = primer_df.drop(rows_to_drop)
    # separate primers into separate dfs of left and right primers
    mask = primer_df['name'].str.contains('LEFT')
    left_primers = primer_df[mask].reset_index(drop=True)
    right_primers = primer_df[~mask].reset_index(drop=True)
    # if alternate primers are included we need to expand the primer dfs so all primers have a mate
    left_primers, right_primers = populate_primer_dfs(left_primers, right_primers)
    if len(left_primers) == len(right_primers):
        return left_primers, right_primers
    else:
        raise ValueError('Left and right primers are different lengths')

def mismatch_positions(left_mismatches,
                    right_mismatches,
                    right_primer_start,
                    left_primer_start):
    """Make primer SNP positions relative to the amplicon"""
    mismatches = {}
    for l_snip in left_mismatches:
        mismatches[int(l_snip + left_primer_start)] = left_mismatches[l_snip]
    for r_snip in right_mismatches:
        mismatches[int(r_snip + right_primer_start)] = right_mismatches[r_snip]
    return mismatches

def reverse_complement(primer_sequence):
    """Get reverse complement of sequence to detect mismatches"""
    base_dict = {'A':'T', 'C':'G', 'G':'C', 'T':'A'}
    complement = ''.join([base_dict[base] for base in primer_sequence])[::-1]
    return complement

def revert_primer(amplicon,
                left_primer,
                left_primer_start,
                left_primer_end,
                right_primer,
                right_primer_start,
                amplicon_primer_mismatches):
    """Swap the primer binding site sequence containing an SNP with \
        with the reference primer sequence"""
    for mismatch in amplicon_primer_mismatches:
        if mismatch - left_primer_start < len(amplicon) / 2:
            reverted = left_primer + amplicon[left_primer_end - left_primer_start:]
        else:
            reverted = amplicon[:right_primer_start - left_primer_start] + reverse_complement(right_primer)
    return reverted

def extract_amplicons(primer_df,
                    sequence_file,
                    random_dropout_probability,
                    match_coverage_mean,
                    match_coverage_sd,
                    mismatch_coverage_mean,
                    mismatch_coverage_sd,
                    dimer_prob,
                    pool1_primers,
                    pool2_primers,
                    seed,
                    output_dir,
                    container_dir):
    """Identify best matching primer pairs and extract the amplicon sequence"""
    # import the simulated squence
    sample_name, sample_sequence = clean_genome(sequence_file)
    # set the seed as the samples name for reproducibility
    rndm = RandomState(int(sample_name.split("_a")[0]))
    # we need to map primers to the simulated genomes so INDELs are properly considered
    aligned_primers, alignment_stats = align_primers(sequence_file,
                                                    primer_df,
                                                    sample_sequence,
                                                    output_dir,
                                                    container_dir)
    # refine primers to the best hits
    left_primers, right_primers = refine_primers(aligned_primers,
                                                alignment_stats)
    # a dictionary to keep track of the frequency and positions of error modes
    amplicon_stats = {}
    # store amplicon coverages for pickling
    amplicon_coverages = {}
    # extract 5'->3' amplicon sequence from sample sequence, including primers
    for position in range(len(left_primers)):
        # assign the amplicon name based on the primer names
        primer_id = left_primers['name'][position] + '---' + right_primers['name'][position]
        # extract the amplicon sequence from the simulated sequence
        amplicon = sample_sequence[left_primers['start'][position]: right_primers['end'][position]]
        total_primer_mismatches = alignment_stats[left_primers['name'][position]]["mismatches"]
        total_primer_mismatches += alignment_stats[right_primers['name'][position]]["mismatches"]
        # get the positions of primer SNPs relative to the amplicon
        amplicon_primer_mismatches = mismatch_positions(alignment_stats[left_primers['name'][position]]["positions"],
                                                        alignment_stats[right_primers['name'][position]]["positions"],
                                                        right_primers['start'][position],
                                                        left_primers['start'][position])
        # populate the amplicon statistics dictionary
        amplicon_stats[primer_id] = {"amplicon_start": int(left_primers['start'][position]),
                                    "left_primer_end": int(left_primers['end'][position]),
                                    "amplicon_end": int(right_primers['end'][position]),
                                    "right_primer_start": int(right_primers['start'][position]),
                                    "has_error": False,
                                    "errors": [],
                                    "primer_mismatches": total_primer_mismatches,
                                    "mismatch_positions": amplicon_primer_mismatches}
        # apply an error mode if there are mismatches unless a SNP occurs in the first or last primer
        if (any(left_primers['name'][position] == name for name in ["SARS-CoV-2_1_LEFT", "nCoV-2019_1_LEFT"]) \
                or any(right_primers['name'][position] == name for name in ["SARS-CoV-2_99_RIGHT", "nCoV-2019_99_RIGHT"])):
            if not total_primer_mismatches == 0:
                new_amplicon = revert_primer(amplicon,
                                            left_primers['seq'][position],
                                            left_primers['start'][position],
                                            left_primers['end'][position],
                                            right_primers['seq'][position],
                                            right_primers['start'][position],
                                            amplicon_primer_mismatches)
                # remove the SNPs in the primer binding site in the truth sequence if SNP is in the first or last primer of the scheme
                with open(sequence_file, "w") as modSeq:
                    modSeq.write(">" + sample_name + "\n" + sample_sequence.replace(amplicon, new_amplicon))
                # continue as if there was no artefact in the amplicon
                amplicon = new_amplicon
                total_primer_mismatches = 0
        else:
            # determine if this amplicon will randomly drop out
            if rndm.binomial(n=1, p=(random_dropout_probability)) == 1:
                amplicon_stats[primer_id]['has_error'] = True
                amplicon_stats[primer_id]["coverage"] = 0
                amplicon_stats[primer_id]["errors"].append("random_dropout")
                amplicon_coverages[primer_id] = 0
                # write out amplicon file
                if not os.path.exists(os.path.join(output_dir, sample_name)):
                    os.mkdir(os.path.join(output_dir, sample_name))
                with open(os.path.join(output_dir, sample_name, primer_id + '.fasta'), 'w') as outFasta:
                    outFasta.write('>' + primer_id + '\n' + amplicon)
                continue
        if not total_primer_mismatches == 0:
            # apply reduced coverage due to SNP or reference primer reversion with 50% probability
            reference_reversion = rndm.binomial(n=1, p=(0.5))
            if not reference_reversion == 1:
                coverage = -1
                while coverage < 0:
                    coverage = int(round(rndm.normal(loc=float(mismatch_coverage_mean),
                                                    scale=float(mismatch_coverage_sd))))
                amplicon_stats[primer_id]["errors"].append("primer_SNP")
            # if the primer reverts to reference, maintain a high sequencing coverage
            else:
                amplicon = revert_primer(amplicon,
                                        left_primers['seq'][position],
                                        left_primers['start'][position],
                                        left_primers['end'][position],
                                        right_primers['seq'][position],
                                        right_primers['start'][position],
                                        amplicon_primer_mismatches)
                coverage = -1
                while coverage < 0:
                    coverage = int(round(rndm.normal(loc=float(match_coverage_mean),
                                                    scale=float(match_coverage_sd))))
                amplicon_stats[primer_id]["errors"].append("primer_reversion")
            amplicon_stats[primer_id]['has_error'] = True
        else:
            # determine sequencing coverage
            coverage = -1
            while coverage < 0:
                coverage = int(round(rndm.normal(loc=float(match_coverage_mean),
                                                scale=float(match_coverage_sd))))
            # check for the possibility of primer dimers and decide if dimer will form with fixed probability
            if left_primers['pool'][position] == "nCoV-2019_1":
                pool_seqs = pool1_primers
            if left_primers['pool'][position] == "nCoV-2019_2":
                pool_seqs = pool2_primers
            dimer = False
            # check if primer dimers can occur between either primers and any others in the pool
            for seq in pool_seqs:
                for dimer_seq in [left_primers["seq"][position], right_primers["seq"][position]]:
                    if dimer_seq[-3:] == reverse_complement(seq[:3]) \
                        and rndm.binomial(n=1, p=(dimer_prob)) == 1:
                        amplicon = dimer_seq[:-3] + reverse_complement(seq)
                        dimer = True
                        if dimer_seq == right_primers["seq"][position]:
                            amplicon = reverse_complement(amplicon)
                    elif seq[-3:] == reverse_complement(dimer_seq[:3]) \
                        and rndm.binomial(n=1, p=(dimer_prob)) == 1:
                        amplicon = seq[-3:] + reverse_complement(dimer_seq)
                        dimer = True
                        if dimer_seq == right_primers["seq"][position]:
                            amplicon = reverse_complement(amplicon)
                    else:
                        pass
            if dimer:
                amplicon_stats[primer_id]['has_error'] = True
                amplicon_stats[primer_id]["errors"].append("primer_dimer")
        amplicon_stats[primer_id]["coverage"] = coverage
        amplicon_coverages[primer_id] = coverage
        # save each amplicon as a separate fasta
        if not os.path.exists(os.path.join(output_dir, sample_name)):
            os.mkdir(os.path.join(output_dir, sample_name))
        with open(os.path.join(output_dir, sample_name, primer_id + '.fasta'), 'w') as outFasta:
            outFasta.write('>' + primer_id + '\n' + amplicon)
    return {str(sample_name): [amplicon_coverages, amplicon_stats]}

def get_options():
    """Parse command line args"""
    parser = argparse.ArgumentParser(description='Split genomic sequences into amplicons and simulate PCR/sequencing error modes.')
    parser.add_argument("--scheme", dest="primer_scheme", required=True,
                        help="artic primer scheme to use",
                        choices=["V3", "V4", "V4.1"],
                        type=str)
    parser.add_argument("--scheme-dir", dest="scheme_dir", required=True,
                        help="artic primer scheme directory",
                        type=str)
    parser.add_argument("--input-dir", dest="input_dir", required=True,
                        help="directory containing simulated sequences",
                        type=str)
    parser.add_argument("--output-dir", dest="output_dir", required=True,
                        help="name of output directory",
                        type=str)
    parser.add_argument("--container-dir", dest="container_dir", required=True,
                        help="name of directory of singularity containers",
                        type=str)
    parser.add_argument("--seed", dest="seed", required=True,
                        help="seed for simulations",
                        type=int)
    parser.add_argument("--dropout-prob", dest="random_dropout_probability", required=True,
                        help="probability of amplicon random dropout",
                        type=float)
    parser.add_argument("--dimer-prob", dest="primer_dimer_probability", required=True,
                        help="probability of primer dimers if pooled primer sequences overlap at the ends",
                        type=float)
    parser.add_argument("--match-mean", dest="match_coverage_mean", required=True,
                        help="mean sequencing coverage for amplicons with exactly matching primers",
                        type=float)
    parser.add_argument("--match-sd", dest="match_coverage_sd", required=True,
                        help="standard deviation of sequencing coverages for amplicons with exactly matching primers",
                        type=float)
    parser.add_argument("--mismatch-mean", dest="mismatch_coverage_mean", required=True,
                        help="mean sequencing coverage for amplicons with mismatching primers",
                        type=float)
    parser.add_argument("--mismatch-sd", dest="mismatch_coverage_sd", required=True,
                        help="standard deviation of sequencing coverages for amplicons with mismatching primers",
                        type=float)
    parser.add_argument("--threads", dest="threads", required=True,
                        help="number of threads to use",
                        type=int)
    args = parser.parse_args()
    return args

def main():
    # parse commmand line args
    args = get_options()
    # import correct primers for primer scheme
    primer_df, pool1_primers, pool2_primers = find_primer_scheme(args.primer_scheme,
                                                                "primer_schemes")
    # create output directory
    if not os.path.exists(args.output_dir):
        os.mkdir(args.output_dir)
    # split input genomic sequences into amplicons by left primer start and right primer end
    sequence_files = glob.glob(os.path.join(args.input_dir, '*.fasta'))
    sys.stderr.write('\nAligning primers and splitting simulated sequences into amplicons\n')
    sample_coverages = {}
    error_stats = {}
    # parallelise amplicon splitting
    job_list = [
        sequence_files[i:i + args.threads] for i in range(0, len(sequence_files), args.threads)
    ]
    for subset in tqdm(job_list):
        amplicon_results = Parallel(n_jobs=args.threads, prefer="threads")(delayed(extract_amplicons)(primer_df,
                                                                                                    genomic_sequence,
                                                                                                    args.random_dropout_probability,
                                                                                                    args.match_coverage_mean,
                                                                                                    args.match_coverage_sd,
                                                                                                    args.mismatch_coverage_mean,
                                                                                                    args.mismatch_coverage_sd,
                                                                                                    args.primer_dimer_probability,
                                                                                                    pool1_primers,
                                                                                                    pool2_primers,
                                                                                                    args.seed,
                                                                                                    args.output_dir,
                                                                                                    args.container_dir) for genomic_sequence in subset)
        for elem in amplicon_results:
            for sample in elem:
                sample_coverages[sample] = elem[sample][0]
                error_stats[sample] = elem[sample][1]
    # Pickle the output
    sys.stderr.write("\nPickling the output\n")
    with open(os.path.join(args.output_dir, 'amplicon_coverages.pickle'), 'wb') as ampOut:
        # Pickle using the highest protocol available.
        pickle.dump(sample_coverages, ampOut, pickle.HIGHEST_PROTOCOL)
    # save amplicon error statistics
    sys.stderr.write("\nPickling the error statistics\n")
    with open(os.path.join(args.output_dir, 'amplicon_statistics.pickle'), 'wb') as statOut:
        pickle.dump(error_stats, statOut, pickle.HIGHEST_PROTOCOL)
    with open(os.path.join(args.output_dir, 'amplicon_statistics.json'), 'w') as jsonOut:
        jsonOut.write(json.dumps(error_stats))
    sys.exit(0)

if __name__== '__main__':
    main()
