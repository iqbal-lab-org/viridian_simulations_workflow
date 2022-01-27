import argparse
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
    parser.add_argument("--genomic-sequence", dest="genomic_sequence", required=True,
                        type=str)
    parser.add_argument("--primer-fasta", dest="primer_fasta", required=True,
                        type=str)
    parser.add_argument("--sai-file", dest="sai_file", required=True,
                        type=str)
    parser.add_argument("--alignment-file", dest="alignment_file", required=True,
                        type=str)
    parser.add_argument("--bam-file", dest="bam_file", required=True,
                        type=str)
    parser.add_argument("--bed-file", dest="bed_file", required=True,
                        type=str)
    args = parser.parse_args()
    return args

def main():
    # parse commmand line args
    args = get_options()
    # aln primer to genomic sequence using bwa, then convert to sam, bam and bed files
    map_command = args.bwa_path + ' index ' + args.genomic_sequence + ' && ' + args.bwa_path + ' aln -n 3 '
    map_command += args.genomic_sequence + ' ' + args.primer_fasta + ' > ' + args.sai_file
    map_command += '; ' + args.bwa_path + ' samse ' + args.genomic_sequence
    map_command += ' ' + args.sai_file +  ' ' + args.primer_fasta
    map_command += ' > ' + args.alignment_file + ' && ' + args.samtools_path + ' view -Sb ' + args.alignment_file + ' > ' + args.bam_file
    map_command += ' && ' + args.bedtools_path + ' -i ' + args.bam_file + ' > ' + args.bed_file
    subprocess.run(map_command, shell=True, check=True)
    sys.exit(0)

if __name__== '__main__':
    main()