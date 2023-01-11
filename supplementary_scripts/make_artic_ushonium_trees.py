import glob
import os
import subprocess
from tqdm import tqdm

def call_ushonium(sampleTsv):
	output_dir = os.path.join(os.path.dirname(sampleTsv), "ushonium_output")
	threads = 64
	# run ushonium
	unshonium_command = 'singularity exec singularity/ushonium/ushonium.img make_jsonl.py --title "ARTIC: ' + os.path.basename(os.path.dirname(output_dir)) + '" --cpus ' + str(threads) + ' ' + sampleTsv + ' ' + output_dir
	subprocess.run(unshonium_command, shell=True, check=True)

# list the tsv
sampleTSVs = glob.glob("artic_ushonium_trees/*/samples.tsv")
for t in tqdm(sampleTSVs):
	call_ushonium(t)
