import argparse
from ctypes import alignment
import os
import subprocess
import sys

def get_options():
    """Parse command line args"""
    parser = argparse.ArgumentParser(description='Align primer sequences to simulated genome and convert to readable formats')
    parser.add_argument("--bwa-path", dest="bwa_path", required=True,
                        type=str)
    parser.add_argument("--samtools-path", dest="samtools_path", required=True,
                        type=str)
    parser.add_argument("--bedtools-path", dest="bedtools_path", required=True,
                        type=str)
    parser.add_argument("--dirs", dest="temp_dirs", required=True,
                        type=str)
    parser.add_argument("--genomic-sequence", dest="genomic_sequence", required=True,
                        type=str)
    args = parser.parse_args()
    return args

def main():
    # parse commmand line args
    args = get_options()
    # import list of amplicons
    with open(os.path.join(args.temp_dirs, "amplicons.txt"), "r") as inDirs:
        amplicon_dirs = inDirs.read().split("\n")
    # aln primer to genomic sequence using bwa, then convert to sam, bam and bed files
    for path in amplicon_dirs:
        primer_fasta = path + ".fasta"
        sai_file = path + ".sai"
        alignment_file = path + ".sam"
        bam_file = path + ".bam"
        bed_file = path + ".bed"
        # align primers to simulated sequence and extract mapped information
        map_command = args.bwa_path + ' index ' + args.genomic_sequence + ' && ' + args.bwa_path + ' aln -n 4 '
        map_command += args.genomic_sequence + ' ' + primer_fasta + ' > ' + sai_file
        map_command += '; ' + args.bwa_path + ' samse ' + args.genomic_sequence
        map_command += ' ' + sai_file +  ' ' + primer_fasta
        map_command += ' > ' + alignment_file + ' && ' + args.samtools_path + ' view -Sb ' + alignment_file + ' > ' + bam_file
        map_command += ' && ' + args.bedtools_path + ' -i ' + bam_file + ' > ' + bed_file
        subprocess.run(map_command, shell=True, check=True)
    sys.exit(0)

if __name__== '__main__':
    main()