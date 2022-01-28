import os
import pandas as pd
import sys
# ensure tests are being run from the top level
sys.path.insert(0, os.getcwd())

from scripts.error_modes import find_primer_scheme, clean_genome, correct_alignment, populate_primer_dfs, refine_primers

# write results for find_primer_scheme
# V3
primer_df, pool1_primers, pool2_primers = find_primer_scheme("V3")
primer_df.to_csv("tests/truth_outputs/V3_scheme.csv", index=False)
with open("tests/truth_outputs/V3_pool1.txt", "w") as o:
    o.write("\n".join(pool1_primers))
with open("tests/truth_outputs/V3_pool2.txt", "w") as o:
    o.write("\n".join(pool2_primers))
# V4
primer_df, pool1_primers, pool2_primers = find_primer_scheme("V4")
primer_df.to_csv("tests/truth_outputs/V4_scheme.csv", index=False)
with open("tests/truth_outputs/V4_pool1.txt", "w") as o:
    o.write("\n".join(pool1_primers))
with open("tests/truth_outputs/V4_pool2.txt", "w") as o:
    o.write("\n".join(pool2_primers))
# V4.1
primer_df, pool1_primers, pool2_primers = find_primer_scheme("V4.1")
primer_df.to_csv("tests/truth_outputs/V4.1_scheme.csv", index=False)
with open("tests/truth_outputs/V4.1_pool1.txt", "w") as o:
    o.write("\n".join(pool1_primers))
with open("tests/truth_outputs/V4.1_pool2.txt", "w") as o:
    o.write("\n".join(pool2_primers))

# write results for correct_alignment
sample_name, sample_sequence =  clean_genome("tests/test_inputs/test.fasta")
bed_content = [None, 16483, 16495]
primer_seq = "TGTGTGCTAATGGACAAGTTTTTGG"
sam_tags, bed_content = correct_alignment(sample_sequence, bed_content, primer_seq)

# write results for populate_primer_dfs
primer_df, pool1_primers, pool2_primers = find_primer_scheme("V4.1")
mask = primer_df['name'].str.contains('LEFT')
left_primers = primer_df[mask].reset_index()
right_primers = primer_df[~mask].reset_index()
left_primers, right_primers = populate_primer_dfs(left_primers, right_primers)
left_primers.to_csv("tests/truth_outputs/left_populated.csv", index=False)
right_primers.to_csv("tests/truth_outputs/right_populated.csv", index=False)

# write results for refine_primers
primer_dict = {"name": ["SARS-CoV-2_10_LEFT", "SARS-CoV-2_10_LEFT_alt1", "SARS-CoV-2_10_RIGHT", "SARS-CoV-2_10_RIGHT_alt1", "SARS-CoV-2_10_RIGHT_alt2"],
               "seq": ["TGAGAAGTGCTCTGCCTATACAGT", "TGAATATCACTTTTGAACTTGATGAAAGGATTG", "TCATCTAACCAATCTTCTTCTTGCTCT", "GGTTGAAGAGCAGCAGAAGTG", "ACCCATAGCAAAAGGTAAAAAGGC"],
               "start": [20, 24, 50, 51, 70],
               "end": [30, 35, 65, 67, 90],
               "mismatches": [3, 1, 1, 2, 0]}
alignment_stats = {"SARS-CoV-2_10_LEFT": {"mismatches": 3}, "SARS-CoV-2_10_LEFT_alt1": {"mismatches": 1}, "SARS-CoV-2_10_RIGHT": {"mismatches": 1}, "SARS-CoV-2_10_RIGHT_alt1": {"mismatches": 2}, "SARS-CoV-2_10_RIGHT_alt2": {"mismatches": 0}}
primer_df = pd.DataFrame(primer_dict)
left_primers, right_primers = refine_primers(primer_df, alignment_stats)
left_primers.to_csv("tests/truth_outputs/left_refined.csv", index=False)
right_primers.to_csv("tests/truth_outputs/right_refined.csv", index=False)

# write results for mismatch_positions
# write results for reverse_complement
# write results for revert_primer
# write results for extract_amplicons
