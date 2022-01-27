import argparse
import glob
from joblib import Parallel, delayed
import json
import numpy as np
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

def find_primer_scheme(scheme):
    """Determine the primer scheme and load the metadata"""
    if scheme == "V3":
        # path of V3 artic primer metadata
        primer_sequences = "V3/nCoV-2019.tsv"
        # import tsv file as pandas dataframe
        primer_df = pd.read_csv(primer_sequences, sep='\t')
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
        primer_sequences = os.path.join(scheme, "SARS-CoV-2.primer.bed")
        # import bed file as string
        with open(primer_sequences, "r") as inV4:
            primers = inV4.read().splitlines()
        # keep track of which primers are in pool1 or pool2
        pool1_primers = []
        pool2_primers = []
        # extract primer names and sequences and convert to df
        primer_dict = {"name": [], "seq": [], "pool": []}
        for row in primers:
            split = row.split("\t")
            primer_dict["name"].append(split[3])
            primer_dict["seq"].append(split[6])
            if split[4] == "1":
                pool1_primers.append(primer_dict["seq"])
                primer_dict["pool"].append("nCoV-2019_1")
            if split[4] == "2":
                pool2_primers.append(primer_dict["seq"])
                primer_dict["pool"].append("nCoV-2019_2")
        if scheme == "V4":
            return pd.DataFrame(primer_dict), pool1_primers, pool2_primers
        elif scheme == "V4.1":
            with open(os.path.join(scheme, "SARS-CoV-2.primer_pairs.tsv")) as inTSV:
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
                    primer_data.append((prospective_name, primer_dict["seq"][primer_dict["name"].index(prospective_name)]))
            primer_df = pd.DataFrame(primer_data, columns=['name','seq'])
            primer_df["pool"] = primer_pools
            return primer_df, pool1_primers, pool2_primers
        else:
            ValueError('Only "V3", "V4" and "V4.1" primer schemes are supported')
    else:
        raise ValueError('Only "V3", "V4" and "V4.1" primer schemes are supported')

def correct_alignment(sample_sequence,
                      bed_content,
                      primer_seq,
                      primer_reversed):
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
                  output_dir):
    """Align primers to genome using BWA and return start and end positions"""
    aligned_starts = []
    aligned_ends = []
    alignment_stats = {}
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
        map_command = 'singularity run singularity/images/BWA.img index ' + genomic_sequence + ' && singularity run singularity/images/BWA.img aln -n 3 '
        map_command += genomic_sequence + ' ' + primer_fasta + ' > ' + os.path.join(temp_dir, row['name'] + '.sai')
        map_command += '; singularity run singularity/images/BWA.img samse ' + genomic_sequence
        map_command += ' ' + os.path.join(temp_dir, row['name'] + '.sai') +  ' ' + primer_fasta
        map_command += ' > ' + alignment_file + ' && singularity run singularity/images/Samtools.img view -Sb ' + alignment_file + ' > ' + bam_file
        map_command += ' && singularity run singularity/images/Bedtools.img bamToBed -i ' + bam_file + ' > ' + bed_file
        subprocess.run(map_command, shell=True, check=True)
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
                                                          mapped.query_sequence,
                                                          mapped.is_reverse)
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
            alignment_stats[row['name']] = {"mismatches": mismatches,
                                            "positions": mismatch_dict}
        # extract start and end position of aligned primer
        aligned_starts.append(int(bed_content[1]))
        aligned_ends.append(int(bed_content[2]))
        # remove the temporary directory
        try:
            shutil.rmtree(temp_dir)
        except FileNotFoundError:
            pass
    # store the mapped stop and start positions
    primer_df['start'] = aligned_starts
    primer_df['end'] = aligned_ends
    return primer_df, alignment_stats

def clean_genome(fasta_file):
    """Split input fasta into sample name and sequence"""
    with open(fasta_file, 'r') as g:
        fasta = g.read().splitlines()
    # split fasta into sample name and sequence
    sample_name = fasta[0].split('>')[1]
    sample_sequence = ''.join(fasta[1:])
    return sample_name, sample_sequence

def refine_primers(primer_df,
                   alignment_stats):
    """This function selects the best matching left and right primer for each amplicon \
        if both match, then the non-alternative primer is selected"""
    rows_to_drop = []
    for name in range(len(primer_df["name"])-1):
        if "alt" in primer_df["name"][name+1] and primer_df["start"][name+1] < primer_df["start"][name] + len(primer_df["seq"][name]):
            if alignment_stats[primer_df["name"][name]]["mismatches"] <= alignment_stats[primer_df["name"][name + 1]]["mismatches"]:
                # drop the alternate primer if it is not the best match
                rows_to_drop.append(name+1)
            else:
                rows_to_drop.append(name)
    primer_df = primer_df.drop(rows_to_drop)
    # separate primers into separate dfs of left and right primers
    mask = primer_df['name'].str.contains('LEFT')
    left_primers = primer_df[mask].reset_index()
    right_primers = primer_df[~mask].reset_index()
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
                      output_dir):
    """Identify best matching primer pairs and extract the amplicon sequence"""
    # set the seed
    np.random.seed(seed)
    # import the simulated squence
    sample_name, sample_sequence = clean_genome(sequence_file)
    # we need to map primers to the simulated genomes so INDELs are properly considered
    aligned_primers, alignment_stats = align_primers(sequence_file,
                                                     primer_df,
                                                     sample_sequence,
                                                     output_dir)
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
                                    "amplicon_end": int(right_primers['end'][position]),
                                    "has_error": False,
                                    "errors": [],
                                    "primer_mismatches": total_primer_mismatches,
                                    "mismatch_positions": amplicon_primer_mismatches}
        # determine if this amplicon will randomly drop out
        if np.random.binomial(n=1, p=(random_dropout_probability)) == 1:
            amplicon_stats[primer_id]['has_error'] = True
            amplicon_stats[primer_id]["coverage"] = 0
            amplicon_stats[primer_id]["errors"].append("random_dropout")
            amplicon_coverages[primer_id] = 0
            continue
        # determine if it's necessary to apply an error mode
        if total_primer_mismatches == 0:
            # determine sequencing coverage
            coverage = -1
            while coverage < 0:
                coverage = int(round(np.random.normal(loc=float(match_coverage_mean),
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
                        and np.random.binomial(n=1, p=(dimer_prob)) == 1:
                        amplicon = dimer_seq[:-3] + reverse_complement(seq)
                        dimer = True
                    elif seq[-3:] == reverse_complement(dimer_seq[:3]) \
                        and np.random.binomial(n=1, p=(dimer_prob)) == 1:
                        amplicon = seq[-3:] + reverse_complement(dimer_seq)
                        dimer = True
                    else:
                        pass
            if dimer:
                amplicon_stats[primer_id]['has_error'] = True
                amplicon_stats[primer_id]["errors"].append("primer_dimer")
        else:
            # apply reduced coverage due to SNP or reference primer reversion with 50% probability
            reference_reversion = np.random.binomial(n=1, p=(0.5))
            if not reference_reversion == 1:
                coverage = -1
                while coverage < 0:
                    coverage = int(round(np.random.normal(loc=float(mismatch_coverage_mean),
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
                    coverage = int(round(np.random.normal(loc=float(match_coverage_mean),
                                                    scale=float(match_coverage_sd))))
                amplicon_stats[primer_id]["errors"].append("primer_reversion")
            amplicon_stats[primer_id]['has_error'] = True
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
    parser.add_argument("--input-dir", dest="input_dir", required=True,
                        help="directory containing simulated sequences",
                        type=str)
    parser.add_argument("--output-dir", dest="output_dir", required=True,
                        help="name of output directory",
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
    primer_df, pool1_primers, pool2_primers = find_primer_scheme(args.primer_scheme)
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
                                                                                                      args.output_dir) for genomic_sequence in subset)
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
