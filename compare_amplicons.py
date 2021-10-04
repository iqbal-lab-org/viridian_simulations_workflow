"""Script to generate visualisations for truth vs viridian-assembled amplicon sequences"""
from collections import Counter
import glob
import matplotlib.pyplot as plt
import os
import sys

def trim_genome(primer_positions,
                reference):
    with open(reference, 'r') as ref:
        ref_seq = ''.join(ref.read().splitlines()[1:])
    scheme_start = int(primer_positions[0].split('\t')[1])
    scheme_end = int(primer_positions[-1].split('\t')[2])
    trimmed_reference = ref_seq[scheme_start: scheme_end]
    return trimmed_reference

def count_errors(truth,
                assembly):
    """Function to find amplicon sequences in viridian assemblies"""
    with open(os.path.join(assembly, 'viridian', 'consensus.final_assembly.fa'), 'r') as con:
        consensus = con.read().replace('>assembly','').replace('\n','')
    if not len(truth) == len(consensus):
        sys.stderr.write('\n[WARNING]: Reference length=' + str(len(truth)) + \
            ', Assembly length=' + str(len(consensus)) + '\n')
    mismatches = 0
    for base in range(len(truth)):
        if not truth[base] == consensus[base]:
            mismatches += 1
    return mismatches

def main():
     # path of V3 artic primer metadata
    primer_positions = "V3/nCoV-2019.primer.bed"
    # import bed file as list
    with open(primer_positions, "r") as pos:
        primer_metadata = pos.read().splitlines()
    # list viridian-generated assemblies
    assemblies = glob.glob('viridian_output' + '/*/')
    # get truth set of amplicons from simulated sequences
    mismatch_counts = []
    for assem in assemblies:
        simulated_genome = os.path.join(input[1], os.path.basename(assem[:-1]), '.fasta')
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

if __name__== '__main__':
    main()