import os
import pandas as pd
import pysam
import subprocess
import sys
# ensure tests are being run from the top level
sys.path.insert(0, os.getcwd())

from scripts.error_modes import find_primer_scheme, clean_genome, correct_alignment, populate_primer_dfs, refine_primers

# uncompress test directories
subprocess.run("tar -zxf tests/test_inputs.tar.gz && tar -zxf tests/truth_outputs.tar.gz", shell=True, check=True)
# test find_primer_scheme
sys.stderr.write("Testing find_primer_scheme\n")
# V3
primer_df, pool1_primers, pool2_primers = find_primer_scheme("V3")
truth_df = pd.read_csv("tests/truth_outputs/V3_scheme.csv")
passed = True
if not primer_df.equals(truth_df):
    passed = False
with open("tests/truth_outputs/V3_pool1.txt", "r") as inPool1:
    truth_pool1 = inPool1.read().split("\n")
if not all(primer in pool1_primers for primer in truth_pool1):
    passed = False
with open("tests/truth_outputs/V3_pool2.txt", "r") as inPool2:
    truth_pool2 = inPool2.read().split("\n")
if not all(primer in pool2_primers for primer in truth_pool2):
    passed = False
# raise error if any of the checks fail
if not passed:
    raise AttributeError("V3 primer scheme test failed")

# V4
primer_df, pool1_primers, pool2_primers = find_primer_scheme("V4")
truth_df = pd.read_csv("tests/truth_outputs/V4_scheme.csv")
passed = True
if not primer_df.equals(truth_df):
    passed = False
with open("tests/truth_outputs/V4_pool1.txt", "r") as inPool1:
    truth_pool1 = inPool1.read().split("\n")
if not all(primer in pool1_primers for primer in truth_pool1):
    passed = False
with open("tests/truth_outputs/V4_pool2.txt", "r") as inPool2:
    truth_pool2 = inPool2.read().split("\n")
if not all(primer in pool2_primers for primer in truth_pool2):
    passed = False
# raise error if any of the checks fail
if not passed:
    raise AttributeError("V4 primer scheme test failed")

# V4.1
primer_df, pool1_primers, pool2_primers = find_primer_scheme("V4.1")
truth_df = pd.read_csv("tests/truth_outputs/V4.1_scheme.csv")
passed = True
if not primer_df.equals(truth_df):
    passed = False
with open("tests/truth_outputs/V4.1_pool1.txt", "r") as inPool1:
    truth_pool1 = inPool1.read().split("\n")
if not all(primer in pool1_primers for primer in truth_pool1):
    passed = False
with open("tests/truth_outputs/V4.1_pool2.txt", "r") as inPool2:
    truth_pool2 = inPool2.read().split("\n")
if not all(primer in pool2_primers for primer in truth_pool2):
    passed = False
# raise error if any of the checks fail
if not passed:
    raise AttributeError("V4.1 primer scheme test failed")
sys.stderr.write("\nTest passed for find_primer_scheme\n")

# test align_primers
sys.stderr.write("\nTesting align_primers\n")
map_command = "singularity run singularity/images/Map.img --dirs tests/test_inputs --genomic-sequence tests/test_inputs/test.fasta"
subprocess.run(map_command, shell=True, check=True)
sys.stderr.write("\nTest passed for align_primers\n")

# test correct_alignment
sys.stderr.write("\nTesting correct_alignment\n")
sample_name, sample_sequence =  clean_genome("tests/test_inputs/test.fasta")
samfile = pysam.AlignmentFile("tests/test_inputs/tmp24g2pmlq/SARS-CoV-2_11_RIGHT.sam", "r")
primer_seq = "GAGGTGTTGCAGGAGCCTTAAA"
sam_tags, bed_content = correct_alignment(sample_sequence, [None, 3470, 3485], primer_seq)
with open("tests/test_inputs/tmp24g2pmlq/SARS-CoV-2_11_RIGHT.bed", 'r') as bedIn:
    true_bed_content = bedIn.read().split('\t')
true_bed_content[0] = None
true_bed_content[1] = int(true_bed_content[1])
true_bed_content[2] = int(true_bed_content[2])
if not sam_tags == [('NM', 0), ('MD', '25')] and not bed_content == true_bed_content[:3]:
    raise AttributeError("correct_alignment test failed")
else:
    sys.stderr.write("\nTest passed for correct_alignment\n")

# test populate_primer_dfs
# separate primers into separate dfs of left and right primers
sys.stderr.write("\nTesting populate_primer_dfs\n")
primer_df, pool1_primers, pool2_primers = find_primer_scheme("V4.1")
mask = primer_df['name'].str.contains('LEFT')
left_primers = primer_df[mask].reset_index()
right_primers = primer_df[~mask].reset_index()
left_primers, right_primers = populate_primer_dfs(left_primers, right_primers)
true_left_primers = pd.read_csv("tests/truth_outputs/left_populated.csv")
true_right_primers = pd.read_csv("tests/truth_outputs/right_populated.csv")
if not left_primers.equals(true_left_primers) or not right_primers.equals(true_right_primers):
    raise AttributeError("populate_primer_dfs test failed")
else:
    sys.stderr.write("\nTest passed for populate_primer_dfs\n")

# test refine_primers
sys.stderr.write("\nTesting refine_primers\n")
primer_dict = {"name": ["SARS-CoV-2_10_LEFT", "SARS-CoV-2_10_LEFT_alt1", "SARS-CoV-2_10_RIGHT", "SARS-CoV-2_10_RIGHT_alt1", "SARS-CoV-2_10_RIGHT_alt2"],
               "seq": ["TGAGAAGTGCTCTGCCTATACAGT", "TGAATATCACTTTTGAACTTGATGAAAGGATTG", "TCATCTAACCAATCTTCTTCTTGCTCT", "GGTTGAAGAGCAGCAGAAGTG", "ACCCATAGCAAAAGGTAAAAAGGC"],
               "start": [20, 24, 50, 51, 70],
               "end": [30, 35, 65, 67, 90],
               "mismatches": [3, 1, 1, 2, 0]}
alignment_stats = {"SARS-CoV-2_10_LEFT": {"mismatches": 3}, "SARS-CoV-2_10_LEFT_alt1": {"mismatches": 1}, "SARS-CoV-2_10_RIGHT": {"mismatches": 1}, "SARS-CoV-2_10_RIGHT_alt1": {"mismatches": 2}, "SARS-CoV-2_10_RIGHT_alt2": {"mismatches": 0}}
primer_df = pd.DataFrame(primer_dict)
left_primers, right_primers = refine_primers(primer_df, alignment_stats)
true_refined_left = pd.read_csv("tests/truth_outputs/left_refined.csv")
true_refined_right = pd.read_csv("tests/truth_outputs/right_refined.csv")
if not left_primers.equals(true_refined_left) or not right_primers.equals(true_refined_right):
    raise AttributeError("refine_primers test failed")
else:
    sys.stderr.write("\nTest passed for refine_primers\n")

# test mismatch_positions
# test reverse_complement
# test revert_primer
# test extract_amplicons

sys.stderr.write("\nAll tests were successful!\n")

# delete outputs and compress test directories
subprocess.run("tar -zcf tests/test_inputs.tar.gz tests/test_inputs/ && tar -zcf tests/truth_outputs.tar.gz tests/truth_outputs/ && rm -rf tests/test_inputs tests/truth_outputs", shell=True, check=True)