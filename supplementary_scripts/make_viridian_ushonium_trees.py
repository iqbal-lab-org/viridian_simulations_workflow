import glob
import os
import subprocess
from tqdm import tqdm

def call_ushonium(sampleTsv):
	output_dir = os.path.join(os.path.dirname(sampleTsv), "ushonium_output")
	threads = 64
	# run ushonium
	unshonium_command = 'singularity run /nfs/research/zi/mhunt/Containers/ushonium.20230303.9efc65b.img --ref_start_end 25,29854 --title "Viridian: ' + os.path.basename(os.path.dirname(output_dir)) + '" --cpus ' + str(threads) + ' ' + sampleTsv + ' ' + output_dir
	subprocess.run(unshonium_command, shell=True, check=True)

# list the tsv
sampleTSVs = glob.glob("ushonium_trees/viridian_ushonium_trees/*/samples.tsv")
for t in tqdm(sampleTSVs):
	if not os.path.exists(os.path.join(os.path.dirname(t), "ushonium_output")):
		call_ushonium(t)
